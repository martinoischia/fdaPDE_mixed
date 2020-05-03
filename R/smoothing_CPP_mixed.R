CPP_smooth.FEM.mixed<-function(locations, observations, num_units, FEMbasis, lambda, covariates = NULL, incidence_matrix = NULL, ndim, mydim, BC = NULL, GCV, GCVMETHOD = 2, nrealizations = 100, DOF=TRUE, DOF_matrix=NULL, search, bary.locations)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(DOF_matrix))
  {
    DOF_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
  }

  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1
  }

  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(num_units) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  DOF_matrix <- as.matrix(DOF_matrix)
  storage.mode(DOF_matrix) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"

  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <-"integer"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"
  storage.mode(search) <- "integer"

  ## Call C++ function
  bigsol <- .Call("regression_Laplace_mixed", locations, observations, num_units, FEMbasis$mesh, FEMbasis$order,
                  mydim, ndim, lambda, covariates, incidence_matrix, BC$BC_indices, BC$BC_values,
                  GCV, GCVMETHOD, nrealizations, DOF, DOF_matrix, search,  bary.locations, PACKAGE = "fdaPDE")
  return(bigsol)
}
