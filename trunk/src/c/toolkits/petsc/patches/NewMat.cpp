/*!\file:  NewMat.cpp
 * \brief create matrix using the Petsc library
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Petsc includes: */
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>

#include "./petscpatches.h"
#include "../../../shared/shared.h"
#include "../../mpi/issmmpi.h"

/*NewMat(int M,int N){{{*/
Mat NewMat(int M,int N,ISSM_MPI_Comm comm){

	/*output:*/
	Mat outmatrix=NULL;

	/*parameters: */
	double sparsity=.001; //default
	int    m,n;
	int    d_nz,o_nz,nnz;

	/*Determine local sizes: */
	m=DetermineLocalSize(M,comm);
	n=DetermineLocalSize(N,comm);

	nnz=(int)((double)M*(double)N*sparsity); //number of non zeros.
	d_nz=(int)((double)nnz/(double)M/2.0); //number of non zeros per row/2
	o_nz=(int)((double)nnz/(double)M/2.0); //number of non zeros per row/2

	#if _PETSC_MAJOR_ == 3 && _PETSC_MINOR_ > 2
	MatCreateAIJ(comm,m,n,M,N,d_nz,NULL,o_nz,NULL,&outmatrix); 
	#else
	MatCreateMPIAIJ(comm,m,n,M,N,d_nz,NULL,o_nz,NULL,&outmatrix); 
	#endif

	return outmatrix;
}
/*}}}*/
/*NewMat(int M,int N,double sparsity,ISSM_MPI_Comm comm){{{*/
Mat NewMat(int M,int N,double sparsity,ISSM_MPI_Comm comm){

	/*output:*/
	Mat outmatrix=NULL;

	/*parameters: */
	int    m,n;
	int    d_nz,o_nz;
	int    nnz;

	/*Determine local sizes: */
	m=DetermineLocalSize(M,comm);
	n=DetermineLocalSize(N,comm);

	nnz=(int)((double)M*(double)N*sparsity); //number of non zeros.
	d_nz=(int)((double)nnz/(double)M/2.0); //number of non zeros per row/2
	o_nz=(int)((double)nnz/(double)M/2.0); //number of non zeros per row/2

	#if _PETSC_MAJOR_ == 3 && _PETSC_MINOR_ > 2
	if(sparsity==1){
		MatCreateDense(comm,m,n,M,N,NULL,&outmatrix); 
	}
	else{
		MatCreateAIJ(comm,m,n,M,N,d_nz,NULL,o_nz,NULL,&outmatrix); 
	}
	#else
	MatCreateMPIAIJ(comm,m,n,M,N,d_nz,NULL,o_nz,NULL,&outmatrix); 
	#endif

	return outmatrix;
}
/*}}}*/
/*NewMat(int M,int N,int connectivity,int numberofdofspernode){{{*/
Mat NewMat(int M,int N,int connectivity,int numberofdofspernode,ISSM_MPI_Comm comm){

	/*output:*/
	Mat outmatrix=NULL;

	/*parameters: */
	int    m,n;
	int    d_nz,o_nz;

	#if _PETSC_MAJOR_ >= 3 
	#if defined(_HAVE_PETSCDEV_) || _PETSC_MINOR_ >=4
	MatType type;
	#else
	const MatType type;
	#endif
	#else
	MatType type;
	#endif

	/*Determine local sizes: */
	m=DetermineLocalSize(M,comm);
	n=DetermineLocalSize(N,comm);

	/*Figure out number of non zeros per row: */
	d_nz=(int)connectivity*numberofdofspernode/2;
	o_nz=(int)connectivity*numberofdofspernode/2;

	MatCreate(comm,&outmatrix);
	MatSetSizes(outmatrix,m,n,M,N);
	MatSetFromOptions(outmatrix);

	/*preallocation  according to type: */
	MatGetType(outmatrix,&type);

	#if _PETSC_MAJOR_ == 2 
	if((strcmp(type,"mpiaij")==0) || (strcmp(type,"aijmumps")==0)){
		MatMPIAIJSetPreallocation(outmatrix,d_nz,NULL,o_nz,NULL);
	}
	#else
	if((strcmp(type,"mpiaij")==0) || (strcmp(type,"mpidense")==0)){
		MatMPIAIJSetPreallocation(outmatrix,d_nz,NULL,o_nz,NULL);
	}
	#endif

	return outmatrix;
}
/*}}}*/
