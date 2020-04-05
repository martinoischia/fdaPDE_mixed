#define R_VERSION_

#include "fdaPDE.h"
#include "regressionData.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "FPCAData.h"
#include "FPCAObject.h"
#include "solverdefinitions.h"
//#include <chrono>      

#include "mixedFEFPCA.h"
#include "mixedFERegression.h"
#include "mixedFEFPCAfactory.h" 
  
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton(InputHandler &regressionData, SEXP Rmesh)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);
	MixedFERegression<InputHandler, Integrator,ORDER, mydim, ndim> regression(mesh, regressionData);
	regression.apply();

	const std::vector<VectorXr>& solution = regression.getSolution();
	const std::vector<Real>& dof = regression.getDOF();
	const MatrixXr & barycenters = regression.getBarycenters();
	const VectorXi & elementIds = regression.getElementIds();


	//Copy result in R memory                
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 2+5+2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution[0].size(), solution.size()));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, solution.size()));
	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < solution.size(); j++)
	{
		for(UInt i = 0; i < solution[0].size(); i++)
			rans[i + solution[0].size()*j] = solution[j][i];
	}

	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt i = 0; i < solution.size(); i++)
	{
		rans1[i] = dof[i];
	}

	//SEND TREE INFORMATION TO R
	SET_VECTOR_ELT(result, 2, Rf_allocVector(INTSXP, 1)); //tree_header information
	int *rans2 = INTEGER(VECTOR_ELT(result, 2));
	rans2[0] = mesh.getTree().gettreeheader().gettreelev();

	SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
	Real *rans3 = REAL(VECTOR_ELT(result, 3));
	for(UInt i = 0; i < ndim*2; i++)
		rans3[i] = mesh.getTree().gettreeheader().domainorig(i);

	SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt i = 0; i < ndim*2; i++)
		rans4[i] = mesh.getTree().gettreeheader().domainscal(i);


	UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
	SET_VECTOR_ELT(result, 5, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
	int *rans5 = INTEGER(VECTOR_ELT(result, 5));
	for(UInt i = 0; i < num_tree_nodes; i++)
			rans5[i] = mesh.getTree().gettreenode(i).getid();

	for(UInt i = 0; i < num_tree_nodes; i++)
			rans5[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

	for(UInt i = 0; i < num_tree_nodes; i++)
			rans5[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

	SET_VECTOR_ELT(result, 6, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
	Real *rans6 = REAL(VECTOR_ELT(result, 6));
	for(UInt j = 0; j < ndim*2; j++)
	{
		for(UInt i = 0; i < num_tree_nodes; i++)
			rans6[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
	}
	
	//SEND BARYCENTER INFORMATION TO R
	SET_VECTOR_ELT(result, 7, Rf_allocVector(INTSXP, elementIds.rows())); //element id of the locations point (vector)
	int *rans7 = INTEGER(VECTOR_ELT(result, 7));
	for(UInt i = 0; i < elementIds.rows(); i++)
		rans7[i] = elementIds(i);

	SET_VECTOR_ELT(result, 8, Rf_allocMatrix(REALSXP, barycenters.rows(), barycenters.cols())); //barycenter information (matrix)
	Real *rans8 = REAL(VECTOR_ELT(result, 8));
	for(UInt j = 0; j < barycenters.cols(); j++)
	{
		for(UInt i = 0; i < barycenters.rows(); i++)
			rans8[i + barycenters.rows()*j] = barycenters(i,j);
	}



	UNPROTECT(1);
	return(result);
}

template<typename Integrator,UInt ORDER, UInt mydim, UInt ndim>
SEXP FPCA_skeleton(FPCAData &fPCAData, SEXP Rmesh, std::string validation)
{

	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh); //, fPCAData.getSearch());
	std::unique_ptr<MixedFEFPCABase<Integrator, ORDER, mydim, ndim>> fpca = MixedFEFPCAfactory<Integrator, ORDER, mydim, ndim>::createFPCAsolver(validation, mesh, fPCAData);
	fpca->apply();

	const std::vector<VectorXr>& loadings = fpca->getLoadingsMat();
	const std::vector<VectorXr>& scores = fpca->getScoresMat();
	const std::vector<Real>& lambdas = fpca->getLambdaPC();
	const std::vector<Real>& variance_explained = fpca->getVarianceExplained();
	const std::vector<Real>& cumsum_percentage = fpca->getCumulativePercentage();
	const std::vector<Real>& var = fpca->getVar();
	const MatrixXr & barycenters = fpca->getBarycenters();
	const VectorXi & elementIds = fpca->getElementIds();

	//Copy result in R memory  
	SEXP result = NILSXP;
	//result = PROTECT(Rf_allocVector(VECSXP, 7)); //**************will be divided to if/ else if there is saveTreeFlag option
	//### why originally 7? Shouldn't it be 6?
	result = PROTECT(Rf_allocVector(VECSXP, 6+5+2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, loadings[0].size(), loadings.size()));
	SET_VECTOR_ELT(result, 1, Rf_allocMatrix(REALSXP, scores[0].size(), scores.size()));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, lambdas.size()));
	SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, variance_explained.size()));
	SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, cumsum_percentage.size()));
	SET_VECTOR_ELT(result, 5, Rf_allocVector(REALSXP, var.size()));
	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < loadings.size(); j++)
	{
		for(UInt i = 0; i < loadings[0].size(); i++)
			rans[i + loadings[0].size()*j] = loadings[j][i];
	}
	
	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt j = 0; j < scores.size(); j++)
	{
		for(UInt i = 0; i < scores[0].size(); i++)
			rans1[i + scores[0].size()*j] = scores[j][i];
	}
	
	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < lambdas.size(); i++)
	{
		rans2[i] = lambdas[i];
	}

	Real *rans3 = REAL(VECTOR_ELT(result, 3));
	for(UInt i = 0; i < variance_explained.size(); i++)
	{
		rans3[i] = variance_explained[i];
	}
	
	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt i = 0; i < cumsum_percentage.size(); i++)
	{
		rans4[i] = cumsum_percentage[i];
	}
	Real *rans5 = REAL(VECTOR_ELT(result, 5));
	for(UInt i = 0; i < var.size(); i++)
	{
		rans5[i] = var[i];
	}


	//TREE INFORMATION
	SET_VECTOR_ELT(result, 6, Rf_allocVector(INTSXP, 7)); //tree_header information
	int *rans6 = INTEGER(VECTOR_ELT(result, 6));
	rans6[0] = mesh.getTree().gettreeheader().gettreeloc();
	rans6[1] = mesh.getTree().gettreeheader().gettreelev();
	rans6[2] = mesh.getTree().gettreeheader().getndimp();
	rans6[3] = mesh.getTree().gettreeheader().getndimt();
	rans6[4] = mesh.getTree().gettreeheader().getnele();
	rans6[5] = mesh.getTree().gettreeheader().getiava();
	rans6[6] = mesh.getTree().gettreeheader().getiend();

	SET_VECTOR_ELT(result, 7, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
	Real *rans7 = REAL(VECTOR_ELT(result, 7));
	for(UInt i = 0; i < ndim*2; i++)
		rans7[i] = mesh.getTree().gettreeheader().domainorig(i);

	SET_VECTOR_ELT(result, 8, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
	Real *rans8 = REAL(VECTOR_ELT(result, 8));
	for(UInt i = 0; i < ndim*2; i++)
		rans8[i] = mesh.getTree().gettreeheader().domainscal(i);


	UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
	SET_VECTOR_ELT(result, 9, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
	int *rans9 = INTEGER(VECTOR_ELT(result, 9));
	for(UInt i = 0; i < num_tree_nodes; i++)
			rans9[i] = mesh.getTree().gettreenode(i).getid();

	for(UInt i = 0; i < num_tree_nodes; i++)
			rans9[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

	for(UInt i = 0; i < num_tree_nodes; i++)
			rans9[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

	SET_VECTOR_ELT(result, 10, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
	Real *rans10 = REAL(VECTOR_ELT(result, 10));
	for(UInt j = 0; j < ndim*2; j++)
	{
		for(UInt i = 0; i < num_tree_nodes; i++)
			rans10[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
	}
	
	//BARYCENTER INFORMATION
	SET_VECTOR_ELT(result, 11, Rf_allocMatrix(REALSXP, barycenters.rows(), barycenters.cols())); //barycenter information (matrix)
	Real *rans11 = REAL(VECTOR_ELT(result, 11));
	for(UInt j = 0; j < barycenters.cols(); j++)
	{
		for(UInt i = 0; i < barycenters.rows(); i++)
			rans11[i + barycenters.rows()*j] = barycenters(i,j);
	}

	SET_VECTOR_ELT(result, 12, Rf_allocVector(INTSXP, elementIds.rows())); //element id of the locations point (vector)
	int *rans12 = INTEGER(VECTOR_ELT(result, 12));
	for(UInt i = 0; i < elementIds.rows(); i++)
		rans12[i] = elementIds(i);

	UNPROTECT(1);

	return(result);
}


template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
SEXP get_integration_points_skeleton(SEXP Rmesh)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);
	FiniteElement<Integrator,ORDER, mydim, ndim> fe;

	SEXP result;
	PROTECT(result=Rf_allocVector(REALSXP, 2*Integrator::NNODES*mesh.num_elements()));
	for(UInt i=0; i<mesh.num_elements(); i++)
	{
		fe.updateElement(mesh.getElement(i));
		for(UInt l = 0;l < Integrator::NNODES; l++)
		{
			Point p = fe.coorQuadPt(l);
			REAL(result)[i*Integrator::NNODES + l] = p[0];
			REAL(result)[mesh.num_elements()*Integrator::NNODES + i*Integrator::NNODES + l] = p[1];
		}
	}

	UNPROTECT(1);
	return(result);
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim, typename A>
SEXP get_FEM_Matrix_skeleton(SEXP Rmesh, EOExpr<A> oper)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	SpMat AMat;
	Assembler::operKernel(oper, mesh, fe, AMat);

	//Copy result in R memory     
	SEXP result;
	result = PROTECT(Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, AMat.nonZeros() , 2));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, AMat.nonZeros()));

	int *rans = INTEGER(VECTOR_ELT(result, 0));
	Real  *rans2 = REAL(VECTOR_ELT(result, 1));
	UInt i = 0;
	for (UInt k=0; k < AMat.outerSize(); ++k)
		{
			for (SpMat::InnerIterator it(AMat,k); it; ++it)
			{
				//std::cout << "(" << it.row() <<","<< it.col() <<","<< it.value() <<")\n"; 
				rans[i] = 1+it.row();
				rans[i + AMat.nonZeros()] = 1+it.col();
				rans2[i] = it.value();
				i++;
			}
		}
	UNPROTECT(1);
	return(result);
}

extern "C" {

//! This function manages the various options for Spatial Regression, Sangalli et al version
/*!
	This function is then called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\param Rsearch an R-integer to decide the search algorithm type (tree or naive or walking search algorithm).
	\return R-vector containg the coefficients of the solution
*/

SEXP regression_Laplace(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim,
					SEXP Rlambda, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
					SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rsearch)
{
    
    //Set input data
	RegressionData regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, Rcovariates, RincidenceMatrix, RBCIndices, 
			RBCValues, DOF, RGCVmethod, Rnrealizations, Rsearch);
	
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

    if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
    	return(regression_skeleton<RegressionData,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh));
    else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
		return(regression_skeleton<RegressionData,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh));
    else if(regressionData.getOrder()==1 && mydim==2 && ndim==3)
		return(regression_skeleton<RegressionData,IntegratorTriangleP2, 1, 2, 3>(regressionData, Rmesh));
   else if(regressionData.getOrder()==2 && mydim==2 && ndim==3)
		return(regression_skeleton<RegressionData,IntegratorTriangleP4, 2, 2, 3>(regressionData, Rmesh));
	else if(regressionData.getOrder()==1 && mydim==3 && ndim==3)
		return(regression_skeleton<RegressionData,IntegratorTetrahedronP2, 1, 3, 3>(regressionData, Rmesh));

    return(NILSXP);
}


/*!
	This function is then called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RK an R-matrix representing the diffusivity matrix of the model
	\param Rbeta an R-vector representing the advection term of the model
	\param Rc an R-double representing the reaction term of the model
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\param Rsearch an R-integer to decide the search algorithm type (tree or naive or walking search algorithm).
	\return R-vector containg the coefficients of the solution
*/

SEXP regression_PDE(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix,
					SEXP RBCIndices, SEXP RBCValues, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rsearch)
{
	RegressionDataElliptic regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, 
		RBCIndices, RBCValues, DOF, RGCVmethod, Rnrealizations, Rsearch);
	
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	if(regressionData.getOrder() == 1 && ndim==2)
		return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 2 && ndim==2)
		return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 1 && ndim==3)
		return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP2, 1, 2, 3>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 2 && ndim==3)
		return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP4, 2, 2, 3>(regressionData, Rmesh));
	return(NILSXP);
}

/*!
	This function is then called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RK an R object representing the diffusivity tensor of the model
	\param Rbeta an R object representing the advection function of the model
	\param Rc an R object representing the reaction function of the model
	\param Ru an R object representing the forcing function of the model
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\param Rsearch an R-integer to decide the search algorithm type (tree or naive or walking search algorithm).
	\return R-vector containg the coefficients of the solution
*/


SEXP regression_PDE_space_varying(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim,
								SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix,
								SEXP RBCIndices, SEXP RBCValues, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rsearch)
{
    //Set data 
	RegressionDataEllipticSpaceVarying regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Ru, Rcovariates, 
		RincidenceMatrix, RBCIndices, RBCValues, DOF,  RGCVmethod, Rnrealizations, Rsearch);
	
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	if(regressionData.getOrder() == 1 && ndim==2)
		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 2 && ndim==2)
		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 1 && ndim==3)
		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, 2, 3>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 2 && ndim==3)
		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, 2, 3>(regressionData, Rmesh));
	return(NILSXP);
}

//! A function required for anysotropic and nonstationary regression (only 2D)
/*! 
    \return points where the PDE space-varying params are evaluated in the R code
*/
SEXP get_integration_points(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
{
	//Declare pointer to access data from C++
	int order = INTEGER(Rorder)[0];
	
	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];
//Not implemented for ndim==3
    if(order == 1 && ndim ==2)
    	return(get_integration_points_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh));
    else if(order == 2 && ndim==2)
    	return(get_integration_points_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh));
    return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging

SEXP get_FEM_mass_matrix(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
{
	int order = INTEGER(Rorder)[0];
	
	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];
	
	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);

    if(order==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, mass));
	if(order==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, mass));
	return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging
SEXP get_FEM_stiff_matrix(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
{
	int order = INTEGER(Rorder)[0];
	
	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];
	
	typedef EOExpr<Stiff> ETMass;   Stiff EStiff;   ETMass stiff(EStiff);

    if(order==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, stiff));
	if(order==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, stiff));
	return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging
SEXP get_FEM_PDE_matrix(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc,
				   SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF,SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rsearch)
{
	RegressionDataElliptic regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, RBCIndices, 
		RBCValues, DOF, RGCVmethod, Rnrealizations, Rsearch);
	
	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	const Real& c = regressionData.getC();
	const Eigen::Matrix<Real,2,2>& K = regressionData.getK();
	const Eigen::Matrix<Real,2,1>& beta = regressionData.getBeta();

    if(regressionData.getOrder()==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	if(regressionData.getOrder()==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging
SEXP get_FEM_PDE_space_varying_matrix(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru,
		   SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF,SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rsearch)
{
	RegressionDataEllipticSpaceVarying regressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Ru, Rcovariates, 
		RincidenceMatrix, RBCIndices, RBCValues, DOF, RGCVmethod, Rnrealizations, Rsearch);
	
	//Get mydim and ndim
	//UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	const Reaction& c = regressionData.getC();
	const Diffusivity& K = regressionData.getK();
	const Advection& beta = regressionData.getBeta();

    if(regressionData.getOrder()==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	if(regressionData.getOrder()==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	return(NILSXP);
}



//! This function manages the various options for SF-PCA
/*!
	This function is than called from R code.
	\param Rdatamatrix an R-matrix containing the datamatrix of the problem.
	\param Rlocations an R-matrix containing the location of the observations.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param RincidenceMatrix an R-matrix representing the incidence matrix defining regions in the model with areal data
	\param Rmydim an R-integer containing the dimension of the problem we are considering.
	\param Rndim an R-integer containing the dimension of the space in which the location are.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RnPC an R-integer specifying the number of principal components to compute.
	\param Rvalidation an R-string containing the method to use for the cross-validation of the penalization term lambda.
	\param RnFolds an R-integer specifying the number of folds to use if K-Fold cross validation method is chosen.		
	\param RGCVmethod an R-integer specifying if the GCV computation has to be exact(if = 1) or stochastic (if = 2).		
	\param Rnrealizations an R-integer specifying the number of realizations to use when computing the GCV stochastically.
	\param Rsearch an R-integer to decide the search algorithm type (tree or naive or walking search algorithm).
	\return R-vector containg the coefficients of the solution
*/
SEXP Smooth_FPCA(SEXP Rlocations, SEXP RbaryLocations, SEXP Rdatamatrix, SEXP Rmesh, SEXP Rorder, SEXP RincidenceMatrix, SEXP Rmydim, SEXP Rndim, 
	SEXP Rlambda, SEXP RnPC, SEXP Rvalidation, SEXP RnFolds, SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rsearch){
	//Set data                
	FPCAData fPCAdata(Rlocations, RbaryLocations, Rdatamatrix, Rorder, RincidenceMatrix, Rlambda, 
						RnPC, RnFolds, RGCVmethod, Rnrealizations, Rsearch);
	UInt mydim=INTEGER(Rmydim)[0]; 
	UInt ndim=INTEGER(Rndim)[0]; 

	std::string validation=CHAR(STRING_ELT(Rvalidation,0));
	
	if(fPCAdata.getOrder() == 1 && mydim==2 && ndim==2)
		return(FPCA_skeleton<IntegratorTriangleP2, 1, 2, 2>(fPCAdata, Rmesh, validation));
	else if(fPCAdata.getOrder() == 2 && mydim==2 && ndim==2)
		return(FPCA_skeleton<IntegratorTriangleP4, 2, 2, 2>(fPCAdata, Rmesh, validation));
	else if(fPCAdata.getOrder() == 1 && mydim==2 && ndim==3)
		return(FPCA_skeleton<IntegratorTriangleP2, 1, 2, 3>(fPCAdata, Rmesh, validation));
	else if(fPCAdata.getOrder() == 2 && mydim==2 && ndim==3)
		return(FPCA_skeleton<IntegratorTriangleP4, 2, 2, 3>(fPCAdata, Rmesh, validation));
	else if(fPCAdata.getOrder() == 1 && mydim==3 && ndim==3)
		return(FPCA_skeleton<IntegratorTetrahedronP2, 1, 3, 3>(fPCAdata, Rmesh, validation));
	return(NILSXP);
	 }      

}