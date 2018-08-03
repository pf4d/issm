/* \file WriteData.c:
 * \brief: general interface for writing data
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./javascriptio.h"
#include "./../../../c/datastructures/datastructures.h"

/*Primitive data types*/
/*FUNCTION WriteData(IssmPDouble** pmatrix, int* pnel, int* matrix, int M,int N){{{*/
void WriteData(IssmPDouble** pmatrix,int* pnel, int* matrix, int M,int N){

	if(pmatrix && matrix){

		/*Copy matrix: */
		IssmPDouble* dmatrix = xNew<IssmPDouble>(M*N); 
		for (int i=0;i<M*N;i++)dmatrix[i]=(IssmPDouble)matrix[i];
		*pmatrix=dmatrix;
		*pnel=M;
	}
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** pmatrix, int* pM, int* pN, , int* matrix, int M,int N){{{*/
void WriteData(IssmPDouble** pmatrix,int* pM, int* pN, int* matrix, int M, int N){

	if(pmatrix && matrix){

		/*Copy matrix: */
		IssmPDouble* dmatrix = xNew<IssmPDouble>(M*N); 
		for (int i=0;i<M*N;i++)dmatrix[i]=(IssmPDouble)matrix[i];
		*pmatrix=dmatrix;
		*pM=M;
		*pN=N;
	}
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** pmatrix, int* pM, IssmPDouble* pN, , int* matrix, int M,int N){{{*/
void WriteData(IssmPDouble** pmatrix,int* pM, int* pN, IssmPDouble* matrix, int M, int N){

	if(pmatrix && matrix){

		/*Copy matrix: */
		IssmPDouble* dmatrix = xNew<IssmPDouble>(M*N); 
		for (int i=0;i<M*N;i++)dmatrix[i]=matrix[i];
		*pmatrix=dmatrix;
		if(pM)*pM=M;
		if(pN)*pN=N;
	}
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** px, int* pnods, double* vector, int M){{{*/
void WriteData(IssmPDouble** px, int* pnods, double* vector, int M){

	if(px && vector){

		IssmPDouble* dx=xNew<IssmPDouble>(M); 
		for(int i=0;i<M;i++)dx[i]=vector[i];
		*px=dx;
		*pnods=M;
	}
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** px, int* pnods, int* vector, int M){{{*/
void WriteData(IssmPDouble** px, int* pnods, int* vector, int M){

	if(px && vector){

		IssmPDouble* dx=xNew<IssmPDouble>(M); 
		for(int i=0;i<M;i++)dx[i]=(IssmPDouble)vector[i];
		*px=dx;
		*pnods=M;
	}
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble* pdouble, IssmPDouble double){{{*/
void WriteData(IssmPDouble* pdouble, IssmPDouble doublein){

	*pdouble=doublein;
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** pdouble, void* nullptr){{{*/
void WriteData(IssmPDouble** pdouble, void*){
	//do nothing
}
/*}}}*/
/*FUNCTION WriteData(char** pstring, char* string){{{*/
void WriteData(char** pstring, char* stringin){

	char* string=xNew<char>(strlen(stringin)+1);
	xMemCpy<char>(string,stringin,strlen(stringin)+1);

	*pstring=string;
}
/*}}}*/
