/*!\file:  VecMerge.cpp
 * \brief merge vector B into A using partitioning vector A(row_partition_vector)=B;
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

void VecMerge(Vec A, Vec B, double* row_partition_vector,int row_partition_size){

	/*Petsc matrix*/
	int lower_row,upper_row,range;
	int* idxm=NULL;
	double* values=NULL;

	/*Vector sizes: */
	int MB;

	VecGetSize(B,&MB);

	/*If the dimension of the partitioning vector is not the same as that of vector B, we have a problem: */
	if ((row_partition_size !=MB) ){
		_error_("Dimensions of partitioning vector incompatible with dimensions of input vector\n");
	}

	/*Get values from vector B and plug them into vector A, using the partitioning vector*/
	VecGetOwnershipRange(B,&lower_row,&upper_row);
	upper_row--;
	range=upper_row-lower_row+1;

	if (range){
		/*This node owns rows of vector B, get them*/
		idxm=xNew<int>(range);
		values=xNew<double>(range);
		for(int i=0;i<range;i++){
			idxm[i]=lower_row+i;
		}
		VecGetValues(B,range,idxm,values);
		/*Now, modify idxm using the partition vector, and plug values into A*/
		for(int i=0;i<range;i++){
			idxm[i]=int(row_partition_vector[lower_row+i])-1; //-1 because partition vector comes from Matlab, where indices start at 1.
		}
		VecSetValues(A,range,idxm,values,INSERT_VALUES);
	}

	/*Assemble vector*/
	VecAssemblyBegin(A);
	VecAssemblyEnd(A);

	/*Free ressources:*/
	xDelete<int>(idxm);
	xDelete<double>(values);
}
