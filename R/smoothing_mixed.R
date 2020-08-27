smooth.FEM.mixed<-function(locations = NULL, observations, FEMbasis, lambda,
                     covariates, random_effect = NULL, PDE_parameters=NULL, incidence_matrix = NULL,
                     BC = NULL, GCV = FALSE, GCVmethod = "Stochastic", TESTFLAG = FALSE, nrealizations = 100, DOF_matrix=NULL, search = "tree", bary.locations = NULL)
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

  space_varying=checkSmoothingParameters_mixed(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda, covariates=covariates, random_effect=random_effect, incidence_matrix=incidence_matrix, 
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


  checkSmoothingParametersSize_mixed(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda, covariates=covariates, random_effect=random_effect, incidence_matrix=incidence_matrix, 
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

  orig_covariates = covariates
  # Start of implementative covariate conversion
  # find number of statistical units, m
  num_units = dim(observations)[2]
  
  ### Converting into big covariate X with both fixed effect and random effect
  m = num_units #num of statistical units
  p = length(random_effect) #num of random-effect coeff
  q = dim(covariates)[2] #num of common-effect coeff
  N = length(observations) #dim(covariates)[1] should be same
  n = dim(observations)[1] # assumes to have same spatial locations across statistical units

  #transform matrix data to vector data
  observations<-as.vector(observations)

  # convert into implementative ver. (length: (q-p) + m*p)
  matrixV = matrix(0, N, p*m)

  if (p<q && p!=0) { #random-effect as subset
    common_cov = as.matrix(covariates[,-random_effect])
    random_cov = as.matrix(covariates[,random_effect])
    for (i in 1:m) {
      matrixV[((i-1)*n + 1):(i*n), ((i-1)*p + 1):(i*p)] = random_cov[((i-1)*n + 1):(i*n),]
    }
    matrixX= cbind(common_cov,matrixV)   
    covariates = matrixX

  } else if (p==q) { #random-effect as full set
    random_cov = as.matrix(covariates)
    for (i in 1:m) {
    matrixV[((i-1)*n + 1):(i*n), ((i-1)*p + 1):(i*p)] = random_cov[((i-1)*n + 1):(i*n),]
    }
    matrixX= matrixV
    covariates = matrixX

  }
  # else if (p==0) { #no random-effect (keep it as it is)
  #   covariates = as.matrix(covariates)
  # }
  ### End of implementative covariate conversion



  # Start of official covariate conversion (to be used in hypothesis testing)
  matrixV= matrix(0, N, p*(m-1))

  if ((p<q && p!=0) || (p==q)) { #random-effect as subset OR random-effect as full set
    random_cov = as.matrix(orig_covariates[,random_effect])
    for (i in 1:(m-1)) {
      matrixV[((i-1)*n + 1):(i*n), ((i-1)*p + 1):(i*p)] = random_cov[((i-1)*n + 1):(i*n),]
    }

    for (i in 1:n) {
      matrixV[((m-1)*n + i), ] = -rep(random_cov[((m-1)*n + i),], (m-1))
    }
    official_covariates= cbind(orig_covariates,matrixV)
  } else if (p==0) { #no random-effect
    official_covariates = as.matrix(covariates)
  }
  ### End of official covariate conversion

  ################## End checking parameters, sizes and conversion #############################

  if(class(FEMbasis$mesh) == 'mesh.2D' & is.null(PDE_parameters)){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.mixed(locations=locations, observations=observations, num_units=num_units, FEMbasis=FEMbasis, lambda=lambda,
                                  covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                  BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations, TESTFLAG=TESTFLAG)

  } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==FALSE){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.mixed(locations=locations, observations=observations, num_units=num_units, FEMbasis=FEMbasis, lambda=lambda,
                                      PDE_parameters = PDE_parameters,
                                      covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                      BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations, TESTFLAG=TESTFLAG)

  } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==TRUE){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.sv.mixed(locations=locations, observations=observations, num_units=num_units, FEMbasis=FEMbasis, lambda=lambda,
                                         PDE_parameters = PDE_parameters,
                                         covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                         BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations, TESTFLAG=TESTFLAG)
  }
  else if(class(FEMbasis$mesh) == 'mesh.2.5D'){

    bigsol = NULL
    print('C++ Code Execution')
    # if(!is.null(locations))
    #   stop("The option locations!=NULL for manifold domains is currently not implemented")
    bigsol = CPP_smooth.manifold.FEM.mixed(locations=locations, observations=observations, num_units=num_units, FEMbasis=FEMbasis, lambda=lambda, 
                                          covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim, 
                                          BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations, TESTFLAG=TESTFLAG)
  
  }else if(class(FEMbasis$mesh) == 'mesh.3D'){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.volume.FEM.mixed(locations=locations, observations=observations, num_units=num_units, FEMbasis=FEMbasis, lambda=lambda, 
                                        covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim, 
                                        BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations, TESTFLAG=TESTFLAG)
  }

  f = bigsol[[1]][1 : (m*nrow(FEMbasis$mesh$nodes)),]
  g = bigsol[[1]][(m*nrow(FEMbasis$mesh$nodes)+1) : (2*m*nrow(FEMbasis$mesh$nodes)),]

  dof = bigsol[[2]]
  GCV_ = bigsol[[3]]
  bestlambda = bigsol[[4]]+1

  
  # Start of coefficient conversion
  matrixX = matrix(data=bigsol[[5]],nrow=ncol(covariates),ncol=length(lambda)) #implementative ver. (length: (q-p) + m*p)

  if (p != 0) { #exists random-effect

    # convert into official coeff (length: q + m*p)
    if (p<q && p!=0) { #random-effect as subset
      betaPrime <- matrixX[1:(q-p),,drop=FALSE]
      b_iPrime <- matrixX[-(1:(q-p)),,drop=FALSE] #split matrixX into 2 matrices
    } else if (p==q) { #random-effect as full set
      b_iPrime <- matrixX
    }
    
    #convert fixed-effect
    beta = matrix(0,nrow=q,ncol=length(lambda))
    indBeta = 1
    indBi = 1
    
    for (i in 1:q) {
      if (!is.element(i,random_effect)) { #beta as it is
        beta[i,] = betaPrime[indBeta,,drop=FALSE] 
        indBeta=indBeta+1
      } else { #convert beta prime to original prime
        temp = numeric(length(lambda))
        for (j in 1:m) {
          temp =temp + b_iPrime[indBi+(j-1)*p,]
        }
        beta[i,]=temp/m
        indBi=indBi+1
      }
    }
  
    #convert random-effect
    b_i = matrix(0,nrow=m*p,ncol=length(lambda))
    indRanEff=1 #this index will be cycled according to random_effect elements
    
    for (i in 1:(m*p)) {
      b_i[i,] = b_iPrime[i,]-beta[random_effect[ifelse(indRanEff!=0,indRanEff, p)],]
      indRanEff = (indRanEff+1)%%p
    }

    #change the name of the rows
    rname=c()
    for (i in 1:m) {
      temp=paste('b_', as.character(i),sep="")
      for (j in 1:p) {
        temp2=paste(temp, as.character(j),sep="")
        rname=c(rname, temp2)
      }
    }

    rownames(b_i) = rname
    # End of coefficient conversion
    
  } else { #if p==0, no random-effect
    beta = matrixX
    b_i = NULL
  }
    

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


  # Make Functional objects
  fit.FEM.mixed  = FEM.mixed(f, num_units, FEMbasis)
  PDEmisfit.FEM.mixed = FEM.mixed(g, num_units, FEMbasis)


  # Save information of Barycenter
  if (is.null(bary.locations)) {
      bary.locations = list(locations=locations, element_ids = bigsol[[11]], barycenters = bigsol[[12]])    
  }
  class(bary.locations) = "bary.locations"

  

  # for hypothesis testing (coeff: before the conversion)
  if(GCV == TRUE && TESTFLAG == TRUE) {
    # for inference of parameters
    psi = bigsol[[13]]
    R0 = bigsol[[14]]
    R1 = bigsol[[15]]
    
    stderr=sqrt(GCV_*(length(observations)-dof)/length(observations))
    # pure_obs_len = length(which(!is.na(observations)))
    # stderr=sqrt(GCV_*(pure_obs_len-dof)/pure_obs_len)
    test.ingredient=test.ingredient(beta = beta[,bestlambda], b_i = b_i[,bestlambda], covariates = official_covariates, 
                                    psi = psi, R0 = R0, R1= R1, stderr=stderr[bestlambda],
                                    nlocs=n, num_units = num_units, bestlambda = bestlambda)
  }

  # Prepare return list
  reslist = NULL
  if(GCV == TRUE) {
    stderr=sqrt(GCV_*(length(observations)-dof)/length(observations))
    # pure_obs_len = length(which(!is.na(observations)))
    # stderr=sqrt(GCV_*(pure_obs_len-dof)/pure_obs_len)
    reslist=list(fit.FEM.mixed = fit.FEM.mixed, PDEmisfit.FEM.mixed = PDEmisfit.FEM.mixed, 
      beta = beta, b_i = b_i, edf = dof, GCV = GCV_, stderr=stderr, bestlambda = bestlambda, bary.locations = bary.locations)

    if (TESTFLAG == TRUE) {
      reslist=list(fit.FEM.mixed = fit.FEM.mixed, PDEmisfit.FEM.mixed = PDEmisfit.FEM.mixed, 
      beta = beta, b_i = b_i, edf = dof, GCV = GCV_, stderr=stderr, bestlambda = bestlambda, test.ingredient= test.ingredient, bary.locations = bary.locations)
    }

  }else {
    reslist=list(fit.FEM.mixed = fit.FEM.mixed, PDEmisfit.FEM.mixed = PDEmisfit.FEM.mixed, beta = beta, b_i = b_i, bary.locations = bary.locations)
  }


  return(reslist)
}


test.ingredient<-function(beta, b_i, covariates, psi, R0, R1, stderr, nlocs, num_units, bestlambda) {
  P = t(R1)%*%solve(R0)%*%R1
  W = covariates #already in the form of official design matrix

  # step 0: need to exclude from b_i
  if (!is.null(b_i)) {
    exclude_ind = seq(num_units, length(b_i), num_units)
    official_b_i = b_i[-exclude_ind]
  }

  # step 1: find matrix S
  Q=diag(nlocs*num_units)-(W%*%(solve(t(W)%*%W))%*%t(W))
  S=psi%*%(solve(t(psi)%*%Q%*%psi + bestlambda*P))%*%t(psi)%*%Q
  trS = sum(diag(S))

  #step2: sigma^2 estimate
  # f_hat=eval.FEM.mixed(mod2$fit.FEM,
                     # locations = loc)[,minGCVind]
  # f_hat=as.matrix(f_hat)
  # z_hat = as.vector(f_hat + W%*%coeff)
  # qq=dim(W)[2]
  # est_sig2 = sum((z-z_hat)^2)/(nlocs*3-qq-trS)
  est_sig2 = stderr^2

  # step3: Var of beta estimates
  varcov = est_sig2*solve(t(W)%*%W) +
    est_sig2*solve(t(W)%*%W)%*%t(W)%*%S%*%t(S)%*%W%*%solve(t(W)%*%W)

  # step4: change the name row and column
  beta_name=c()
  for (i in 1:length(beta)) {
    temp=paste('beta_', as.character(i),sep="")
    beta_name=c(beta_name, temp)
  }
  names(beta) = beta_name

  if (!is.null(b_i)) {
    varcov_name = c(names(beta), names(official_b_i))
    rownames(varcov) = varcov_name
    colnames(varcov) = varcov_name
    test.ingredient = list(coeff = c(beta, official_b_i), varcov = varcov) 
  } else {
    varcov_name = names(beta)
    rownames(varcov) = varcov_name
    colnames(varcov) = varcov_name
    test.ingredient = list(coeff = beta, varcov = varcov) 
  }

     
  return (test.ingredient)
}
