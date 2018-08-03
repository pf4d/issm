/*\file FetchData.cpp:
 * \brief: general I/O interface to fetch data in javascript
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./javascriptio.h"
#include <cstring> 

/*Primitive data types*/
/*FUNCTION FetchData(char** pstring, char* string){{{*/
void FetchData(char** pstring, char* stringin){

	char* string=xNew<char>(strlen(stringin)+1);
	xMemCpy<char>(string,stringin,strlen(stringin)+1);

	*pstring=string;
}
/*}}}*/
/*FUNCTION FetchData(int* pinteger,int integer){{{*/
void FetchData(int* pinteger,int integer){
	
	*pinteger = integer;
}
/*}}}*/
/*FUNCTION FetchData(double* pscalar,double scalar){{{*/
void FetchData(double* pscalar,double scalar){
	
	*pscalar = scalar;
}
/*}}}*/
/*FUNCTION FetchData(double **pvector, double* vectorin, int nods){{{*/
void FetchData(double** pvector, double* vectorin, int nods){

	double* vector=NULL;
	
	vector=xNew<IssmPDouble>(nods); xMemCpy<IssmPDouble>(vector,vectorin,nods);
	
	*pvector=vector;
}
/*}}}*/
/*FUNCTION FetchData(double **pvector, int* pnods, double* vectorin, int nods){{{*/
void FetchData(double** pvector, int* pnods, double* vectorin, int nods){

	double* vector=NULL;
	
	vector=xNew<IssmPDouble>(nods); xMemCpy<IssmPDouble>(vector,vectorin,nods);
	
	*pvector=vector;
	*pnods=nods;
}
/*}}}*/
/*FUNCTION FetchData(double **pmatrix, int* pM, int* pN, int* matrix, int M, int N){{{*/
void FetchData(double **pmatrix, int* pM, int* pN, int* matrixin, int M, int N){

	double* matrix=NULL;
	
	if(pmatrix && matrixin){ 

		matrix=xNew<IssmPDouble>(M*N); 
		for(int i=0;i<M*N;i++)matrix[i]=(IssmPDouble)matrixin[i];
		if (pM)*pM=M;
		if (pN)*pN=N;
		*pmatrix=matrix;
	}
}
/*}}}*/
/*FUNCTION FetchData(double **pmatrix, int* pM, int* pN, double* matrix, int M, int N){{{*/
void FetchData(double **pmatrix, int* pM, int* pN, double* matrixin, int M, int N){

	double* matrix=NULL;
	
	if(pmatrix && matrixin){ 

		matrix=xNew<IssmPDouble>(M*N); 
		for(int i=0;i<M*N;i++)matrix[i]=matrixin[i];
		if (pM)*pM=M;
		if (pN)*pN=N;
		*pmatrix=matrix;
	}
}
/*}}}*/
/*FUNCTION FetchData(int **pmatrix, int* pM, int* pN, int* matrix, int M, int N){{{*/
void FetchData(int **pmatrix, int* pM, int* pN, int* matrixin, int M, int N){

	int* matrix=NULL;
	
	if(pmatrix && matrixin){ 

		matrix=xNew<int>(M*N);xMemCpy<int>(matrix,matrixin,M*N); 
		if (pM)*pM=M;
		if (pN)*pN=N;
		*pmatrix=matrix;
	}
}
/*}}}*/
/*ISSM objects*/
/*FUNCTION FetchData(Contours** pcontours,double* x, double* y, int nods){{{*/
void FetchData(Contours** pcontours,double* x, double* y, int nods){

	int             numcontours,index,test1,test2;
	char            *contourname = NULL;
	Contours        *contours    = NULL;
	Contour<double> *contouri    = NULL;

	/*only 1 contour for now: */
	contours=new Contours();

	if (nods){
			
		contouri=new Contour<double>();
		contouri->nods=nods;
		contouri->x=xNew<IssmPDouble>(nods); xMemCpy<IssmPDouble>(contouri->x,x,nods);
		contouri->y=xNew<IssmPDouble>(nods); xMemCpy<IssmPDouble>(contouri->y,y,nods);

		contours->AddObject(contouri);
	}
	
	*pcontours=contours;
}
/*}}}*/
/*FUNCTION FetchData(Options** poptions,int NRHS, int nrhs, const char* optionname, double optionvalue){{{*/
void FetchData(Options** poptions,int NRHS, int nrhs, const char* optionname, double optionvalue){

	/*Initialize output*/
	Options* options=new Options();
	
	GenericOption<double> *odouble = NULL;

	/*check and parse the name  */
	odouble=new GenericOption<double>();
	odouble->name =xNew<char>(strlen(optionname)+1);
	memcpy(odouble->name,optionname,(strlen(optionname)+1)*sizeof(char));
	odouble->value=optionvalue;
	odouble->numel=1;
	odouble->ndims=1;
	odouble->size=NULL;
	 
	options->AddOption((Option*)odouble);

	/*Assign output pointers:*/
	*poptions=options;
}
/*}}}*/
