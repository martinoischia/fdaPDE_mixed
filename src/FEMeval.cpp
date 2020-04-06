/*
 * FEMeval.cpp
 *
 *  Created on: Aug 16, 2015
 *      Author: eardi
 */


#define R_VERSION_

#include "fdaPDE.h"
//#include "IO_handler.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "evaluator.h"
#include "projection.h"
#include <chrono>  // for high_resolution_clock

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP tree_mesh_skeleton(SEXP Rmesh) {  
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	//Copy result in R memory                
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 5));
	

	//SEND TREE INFORMATION TO R
	SET_VECTOR_ELT(result, 0, Rf_allocVector(INTSXP, 1)); //tree_header information
	int *rans = INTEGER(VECTOR_ELT(result, 0));
	rans[0] = mesh.getTree().gettreeheader().gettreelev();

	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt i = 0; i < ndim*2; i++)
		rans1[i] = mesh.getTree().gettreeheader().domainorig(i);

	SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < ndim*2; i++)
		rans2[i] = mesh.getTree().gettreeheader().domainscal(i);


	UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
	SET_VECTOR_ELT(result, 3, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
	int *rans3 = INTEGER(VECTOR_ELT(result, 3));
	for(UInt i = 0; i < num_tree_nodes; i++)
			rans3[i] = mesh.getTree().gettreenode(i).getid();

	for(UInt i = 0; i < num_tree_nodes; i++)
			rans3[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

	for(UInt i = 0; i < num_tree_nodes; i++)
			rans3[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

	SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt j = 0; j < ndim*2; j++)
	{
		for(UInt i = 0; i < num_tree_nodes; i++)
			rans4[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
	}
	

	UNPROTECT(1);
	return(result);
}

extern "C" {
//! This function manages the various option for the solution evaluation.
/*!
	This function is then called from R code.
	Calls the walking algoritm for efficient point location inside the mesh in 2D.

	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rlocations an R-matrix (seen as an array) containing the xyz coordinates of the points where the solution has to be evaluated
	\param RincidenceMatrix an R-matrix for the incidence matrix defining the regions in the case of areal data
	\param Rcoef an R-vector the coefficients of the solution
	\param Rorder an R integer containg the order of the solution
	\param Rfast an R integer 0 for Naive location algorithm, 1 for Walking Algorithm (can miss location for non convex meshes)
*/


SEXP eval_FEM_fd(SEXP Rmesh, SEXP Rlocations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rorder, SEXP Rfast, SEXP Rmydim, SEXP Rndim, SEXP Rsearch, SEXP RbaryLocations)
{
	int n_X = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	int nRegions = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
	int nElements = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[1]; //number of triangles/tetrahedron if areal data
	
	std::vector<UInt> element_id;	
	Real **barycenters;

	//RECIEVE BARYCENTER INFORMATION FROM R
	if (TYPEOF(RbaryLocations) != 0) { //have location information
		element_id.assign(INTEGER(VECTOR_ELT(RbaryLocations, 1)), INTEGER(VECTOR_ELT(RbaryLocations, 1))+n_X);
		UInt n_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[0]; //barycenter rows (number of locations)
		UInt p_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[1]; //barycenter columns (number of vertices)

		barycenters = (Real**) malloc(sizeof(Real*)*n_);
		for (int i=0; i<n_; i++)
		{
			barycenters[i] = (Real*) malloc(sizeof(Real)*p_);
			for (int j=0; j<p_; j++)
			{
				barycenters[i][j] = REAL(VECTOR_ELT(RbaryLocations, 2))[i+n_*j];
			}
		}
	}


	//Declare pointer to access data from C++
	double *X, *Y, *Z;
	UInt **incidenceMatrix;
	double *coef;
	int order, mydim, ndim, search;
	bool fast;

	coef  = REAL(Rcoef);
    order = INTEGER(Rorder)[0];
    fast  = INTEGER(Rfast)[0];
    mydim = INTEGER(Rmydim)[0];
    ndim  = INTEGER(Rndim)[0];
	search  = INTEGER(Rsearch)[0];

	X = (double*) malloc(sizeof(double)*n_X);
	Y = (double*) malloc(sizeof(double)*n_X);
	Z = (double*) malloc(sizeof(double)*n_X);
	incidenceMatrix = (UInt**) malloc(sizeof(UInt*)*nRegions);

    // Cast all computation parameters
	if (ndim==3)
	{
		for (int i=0; i<n_X; i++)
		{
			X[i] = REAL(Rlocations)[i + n_X*0];
			Y[i] = REAL(Rlocations)[i + n_X*1];
			Z[i] = REAL(Rlocations)[i + n_X*2];
		}
	}
	else //ndim==2
	{
		for (int i=0; i<n_X; i++)
		{
			X[i] = REAL(Rlocations)[i + n_X*0];
			Y[i] = REAL(Rlocations)[i + n_X*1];
			Z[i] = 0;
		}
	}
	for (int i=0; i<nRegions; i++)
	{
		incidenceMatrix[i] = (UInt*) malloc(sizeof(UInt)*nElements);
		for (int j=0; j<nElements; j++)
		{
			incidenceMatrix[i][j] = INTEGER(RincidenceMatrix)[i+nRegions*j];
		}
	}

    SEXP result;

	if (n_X>0) //pointwise data
	{
		PROTECT(result = Rf_allocVector(REALSXP, n_X));
		std::vector<bool> isinside(n_X);
		//Set the mesh
		if(order==1 && mydim==2 && ndim==2)
		{
			MeshHandler<1,2,2> mesh(Rmesh, search);
			Evaluator<1,2,2> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==2 && mydim==2 && ndim==2)
		{
			MeshHandler<2,2,2> mesh(Rmesh, search);
			Evaluator<2,2,2> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information	
				evaluator.eval(X, Y, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==1 && mydim==2 && ndim==3)
		{ 
			MeshHandler<1,2,3> mesh(Rmesh, search);
			Evaluator<1,2,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}

		}
		else if(order==2 && mydim==2 && ndim==3)
		{
			MeshHandler<2,2,3> mesh(Rmesh, search);
			Evaluator<2,2,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}
		else if(order==1 && mydim==3 && ndim==3)
		{ 
			MeshHandler<1,3,3> mesh(Rmesh, search);
			Evaluator<1,3,3> evaluator(mesh);
			if (TYPEOF(RbaryLocations) == 0) { //doesn't have location information	
				evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
			} else { //have location information
				evaluator.evalWithInfo(X, Y, Z, n_X, coef, fast, REAL(result), isinside, element_id, barycenters);
			}
		}

		for (int i=0; i<n_X; ++i)
		{
			if(!(isinside[i]))
			{
				REAL(result)[i]=NA_REAL;
			}
		}
	}
	else //areal data
	{
		PROTECT(result = Rf_allocVector(REALSXP, nRegions));
		if(order==1 && mydim==2 && ndim==2)
		{
			MeshHandler<1,2,2> mesh(Rmesh);
			Evaluator<1,2,2> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==2 && mydim==2 && ndim==2)
		{
			MeshHandler<2,2,2> mesh(Rmesh);
			Evaluator<2,2,2> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==1 && mydim==2 && ndim==3)
		{ 
			MeshHandler<1,2,3> mesh(Rmesh);
			Evaluator<1,2,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==2 && mydim==2 && ndim==3)
		{
			MeshHandler<2,2,3> mesh(Rmesh);
			Evaluator<2,2,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==1 && mydim==3 && ndim==3)
		{ 
			MeshHandler<1,3,3> mesh(Rmesh);
			Evaluator<1,3,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));

		}
	}

	free(X); free(Y); free(Z);
	for (int i=0; i<nRegions; i++)
	{
		free(incidenceMatrix[i]);
	}
	free(incidenceMatrix);
	

	UNPROTECT(1);
    // result list
    return(result);
}


SEXP points_projection(SEXP Rmesh, SEXP Rlocations)
{
	int n_X = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	//Declare pointer to access data from C++
	double X, Y, Z;

    // Cast all computation parameters
    std::vector<Point> deData_(n_X); // the points to be projected
    std::vector<Point> prjData_(n_X); // the projected points
    
    //RECIEVE PROJECTION INFORMATION FROM R
    for (int i=0; i<n_X; i++)
	{
		X = REAL(Rlocations)[i + n_X*0];
		Y = REAL(Rlocations)[i + n_X*1];
		Z = REAL(Rlocations)[i + n_X*2];
		deData_[i]=Point(X,Y,Z);
	}

    SEXP result;

	if (n_X>0) //pointwise data
	{
		PROTECT(result = Rf_allocMatrix(REALSXP, n_X, 3));
		UInt order = INTEGER(VECTOR_ELT(Rmesh,4))[0];

		if (order == 1) {
		MeshHandler<1,2,3> mesh(Rmesh);
		projection<1,2,3> projector(mesh, deData_);
		prjData_ = projector.computeProjection();
		}

		// if (order == 2) {
		// MeshHandler<2,2,3> mesh(Rmesh);
		// projection<2,2,3> projector(mesh, deData_);
		// prjData_ = projector.computeProjection();
		// }
	}

	for (int i=0; i<n_X; ++i)
	{
		REAL(result)[i + n_X*0]=prjData_[i][0];
		REAL(result)[i + n_X*1]=prjData_[i][1];
		REAL(result)[i + n_X*2]=prjData_[i][2];
	}

	UNPROTECT(1);
    // result matrix
    return(result);
}

SEXP tree_mesh_construction(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim) {
	UInt ORDER=INTEGER(Rorder)[0];
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	if(ORDER == 1 && mydim==2 && ndim==2)
		return(tree_mesh_skeleton<1, 2, 2>(Rmesh));
	else if(ORDER == 2 && mydim==2 && ndim==2)
		return(tree_mesh_skeleton<2, 2, 2>(Rmesh));
	else if(ORDER == 1 && mydim==2 && ndim==3)
		return(tree_mesh_skeleton<1, 2, 3>(Rmesh));
	else if(ORDER == 2 && mydim==2 && ndim==3)
		return(tree_mesh_skeleton<2, 2, 3>(Rmesh));
	else if(ORDER == 1 && mydim==3 && ndim==3)
		return(tree_mesh_skeleton<1, 3, 3>(Rmesh));
	return(NILSXP);
 }

}

