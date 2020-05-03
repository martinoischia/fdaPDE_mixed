smooth.FEM.mixed<-function(locations = NULL, observations, FEMbasis, lambda,
                     covariates = NULL, random_effect = NULL, PDE_parameters=NULL, incidence_matrix = NULL,
                     BC = NULL, GCV = FALSE, GCVmethod = "Stochastic", nrealizations = 100, DOF_matrix=NULL, search = "tree", bary.locations = NULL)
{
  if(class(FEMbasis$mesh) == "mesh.2D"){
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.2.5D"){
    ndim = 3
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.3D"){
    ndim = 3
    mydim = 3
  }else{
    stop('Unknown mesh class')
  }
  ##################### Checking parameters, sizes and conversion ################################

  if(GCVmethod=="Stochastic")
    GCVMETHOD=2
  else if(GCVmethod=="Exact")
    GCVMETHOD=1
  else{
    stop("GCVmethod must be either Stochastic or Exact")
  }

  if(search=="naive")
    search=1
  else if(search=="tree")
    search=2
  else{
    stop("search must be either tree or naive.")
  }

  #if locations is null but bary.locations is not null, use the locations in bary.locations
  if(is.null(locations) & !is.null(bary.locations)) {
    locations = bary.locations$locations
    locations = as.matrix(locations)
  }

  DOF=TRUE
  if(!is.null(DOF_matrix))
    DOF=FALSE

  space_varying=checkSmoothingParameters_mixed(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda, covariates=covariates, incidence_matrix=incidence_matrix, 
    BC=BC, GCV=GCV, PDE_parameters=PDE_parameters, GCVmethod=GCVMETHOD , nrealizations=nrealizations, search=search, bary.locations=bary.locations)
  
  ## Converting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
    observations = as.matrix(observations)
    lambda = as.matrix(lambda)
  if(!is.null(covariates))
    covariates = as.matrix(covariates)
  if(!is.null(DOF_matrix))
    DOF_matrix = as.matrix(DOF_matrix)
  if(!is.null(incidence_matrix))
    incidence_matrix = as.matrix(incidence_matrix)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }

  
  # if I have PDE non-sv case I need (constant) matrices as parameters

  if(!is.null(PDE_parameters) & space_varying==FALSE)
  {
    PDE_parameters$K = as.matrix(PDE_parameters$K)
    PDE_parameters$b = as.matrix(PDE_parameters$b)
    PDE_parameters$c = as.matrix(PDE_parameters$c)
  }


  checkSmoothingParametersSize_mixed(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda, covariates=covariates, incidence_matrix=incidence_matrix, 
    BC=BC, GCV=GCV, space_varying=space_varying, PDE_parameters=PDE_parameters, ndim=ndim, mydim=mydim)
  
  # Check whether the locations coincide with the mesh nodes (should be put after all the validations)
  if (!is.null(locations)) {
    if(dim(locations)[1]==dim(FEMbasis$mesh$nodes)[1] & dim(locations)[2]==dim(FEMbasis$mesh$nodes)[2]) {
      sum1=0
      sum2=0
      for (i in 1:nrow(locations)) {
      sum1 = sum1 + abs(locations[i,1]-FEMbasis$mesh$nodes[i,1])
      sum2 = sum2 + abs(locations[i,2]-FEMbasis$mesh$nodes[i,2])
      }
      if (sum1==0 & sum2==0) {
        message("No search algorithm is used because the locations coincide with the nodes.")
        locations = NULL #In principle, R uses pass-by-value semantics in its function calls. So put ouside of checkSmoothingParameters function.
      }
    }
  }

  print('********Start of covariate conversion')
  # find number of statistical units, m
  num_units = dim(observations)[2]
  
  ### Converting into big covariate X with both fixed effect and random effect
  m = num_units
  p = length(random_effect)
  N = length(observations) #dim(covariates)[1] should be same
  n = dim(observations)[1] # assumes to have same spatial locations across statistical units

  #transform matrix data to vector data
  observations<-as.vector(observations)

  common_cov = as.matrix(covariates[,-random_effect])
  random_cov = as.matrix(covariates[,random_effect])
  
  matrixV = matrix(0, N, p*m)
  
  for (i in 1:m) {
    matrixV[((i-1)*n+1):(i*n), ((i-1)*p+1):(i*p)] = random_cov[((i-1)*n+1):(i*n),]
  }

  matrixX= cbind(common_cov,matrixV)
  covariates = matrixX
  ### End of conversion

  print('********End of covariate conversion')
  ################## End checking parameters, sizes and conversion #############################

  if(class(FEMbasis$mesh) == 'mesh.2D' & is.null(PDE_parameters)){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.mixed(locations=locations, observations=observations, num_units=num_units, FEMbasis=FEMbasis, lambda=lambda,
                                  covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                  BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations)

    numnodes = nrow(FEMbasis$mesh$nodes)

  } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==FALSE){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.mixed(locations=locations, observations=observations, num_units=num_units, FEMbasis=FEMbasis, lambda=lambda,
                                      PDE_parameters = PDE_parameters,
                                      covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                      BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations)

    numnodes = nrow(FEMbasis$mesh$nodes)

  } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==TRUE){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.sv.mixed(locations=locations, observations=observations, num_units=num_units, FEMbasis=FEMbasis, lambda=lambda,
                                         PDE_parameters = PDE_parameters,
                                         covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                         BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations)

    numnodes = nrow(FEMbasis$mesh$nodes)

  }
  # else if(class(FEMbasis$mesh) == 'mesh.2.5D'){

  #   bigsol = NULL
  #   print('C++ Code Execution')
  #   # if(!is.null(locations))
  #   #   stop("The option locations!=NULL for manifold domains is currently not implemented")
  #   bigsol = CPP_smooth.manifold.FEM.mixed(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda, 
  #                                         covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim, 
  #                                         BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations)

  #   numnodes = FEMbasis$mesh$nnodes

  # }else if(class(FEMbasis$mesh) == 'mesh.3D'){

  #   bigsol = NULL
  #   print('C++ Code Execution')
  #   bigsol = CPP_smooth.volume.FEM.mixed(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda, 
  #                                       covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim, 
  #                                       BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations)

  #   numnodes = FEMbasis$mesh$nnodes
  # }

  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]

  dof = bigsol[[2]]
  GCV_ = bigsol[[3]]
  bestlambda = bigsol[[4]]+1

  if(!is.null(covariates))
    beta = matrix(data=bigsol[[5]],nrow=ncol(covariates),ncol=length(lambda))
  else
    beta = NULL

   # Save information of Tree Mesh
    tree_mesh = list(
    treelev = bigsol[[6]][1],
    header_orig= bigsol[[7]], 
    header_scale = bigsol[[8]],
    node_id = bigsol[[9]][,1],
    node_left_child = bigsol[[9]][,2],
    node_right_child = bigsol[[9]][,3],
    node_box= bigsol[[10]])

  # Reconstruct FEMbasis with tree mesh
  mesh.class= class(FEMbasis$mesh)
  if (is.null(FEMbasis$mesh$treelev)) { #if doesn't exist the tree information
    FEMbasis$mesh = append(FEMbasis$mesh, tree_mesh)
  } #if already exist the tree information, don't append
  class(FEMbasis$mesh) = mesh.class  

  # Save information of Barycenter
  if (is.null(bary.locations)) {
      bary.locations = list(locations=locations, element_ids = bigsol[[11]], barycenters = bigsol[[12]])    
  }
  class(bary.locations) = "bary.locations"

  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)

  # Prepare return list
  reslist = NULL

  if(GCV == TRUE)
  {
    stderr=sqrt(GCV_*(length(observations)-dof)/length(observations))
    reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM,
            beta = beta, edf = dof, GCV = GCV_, stderr=stderr, bestlambda = bestlambda, bary.locations = bary.locations)
  }else{
    reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta = beta, bary.locations = bary.locations)
  }

  return(reslist)
}