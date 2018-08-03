/*\file FetchData.cpp:
 * \brief: general I/O interface to fetch data in matlab
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./matlabio.h"
#include <cstring> 

/*Primitive data types*/
/*FUNCTION FetchData(double** pmatrix,int* pM,int *pN,const mxArray* dataref){{{*/
void FetchData(double** pmatrix,int* pM,int *pN,const mxArray* dataref){

	double*  outmatrix=NULL;
	int      outmatrix_rows,outmatrix_cols;

	if(mxIsEmpty(dataref) ){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix_cols=0;
		outmatrix=NULL;
	}
	else if( mxIsClass(dataref,"double") || 
				mxIsClass(dataref,"single") || 
				mxIsClass(dataref,"int16") || 
				mxIsClass(dataref,"int8") || 
				mxIsClass(dataref,"uint8")){
		/*Check dataref is not pointing to NaN: */
		if ( mxIsNaN(*(mxGetPr(dataref))) && (mxGetM(dataref)==1) && (mxGetN(dataref)==1) ){
			outmatrix_rows=0;
			outmatrix_cols=0;
			outmatrix=NULL;
		}
		else{
			if(!mxIsClass(dataref,"double") && !mxIsClass(dataref,"single")){
				_printf_("Warning: converting matlab data from '" << mxGetClassName(dataref) << "' to 'double'\n");
			}
			/*Convert matlab matrix to double* matrix: */
			MatlabMatrixToDoubleMatrix(&outmatrix,&outmatrix_rows,&outmatrix_cols,dataref);
		}
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM)*pM=outmatrix_rows;
	if (pN)*pN=outmatrix_cols;

}
/*}}}*/
/*FUNCTION FetchData(double** pmatrix,int* pnumel,int* pndims,int** psize,const mxArray* dataref){{{*/
void FetchData(double** pmatrix,int* pnumel,int* pndims,int** psize,const mxArray* dataref){

	int     outmatrix_numel,outmatrix_ndims;
	double *outmatrix       = NULL;
	int    *outmatrix_size  = NULL;

	if(mxIsEmpty(dataref) ){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_numel = 0;
		outmatrix_ndims = 0;
		outmatrix_size  = NULL;
		outmatrix       = NULL;
	}
	else if( mxIsClass(dataref,"double") || 
				mxIsClass(dataref,"single") || 
				mxIsClass(dataref,"int16") || 
				mxIsClass(dataref,"int8") || 
				mxIsClass(dataref,"uint8")){

		/*Check dataref is not pointing to NaN: */
		if (mxIsNaN(*(mxGetPr(dataref))) && (mxGetNumberOfElements(dataref)==1)){
			outmatrix_numel = 0;
			outmatrix_ndims = 0;
			outmatrix_size  = NULL;
			outmatrix       = NULL;
		}
		else{
			if(!mxIsClass(dataref,"double") && !mxIsClass(dataref,"single")){
				_printf_("Warning: converting matlab data from '" << mxGetClassName(dataref) << "' to 'double'\n");
			}
			/*Convert matlab n-dim array to double* matrix: */
			MatlabNArrayToNArray(&outmatrix,&outmatrix_numel,&outmatrix_ndims,&outmatrix_size,dataref);
		}
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pnumel)*pnumel=outmatrix_numel;
	if (pndims)*pndims=outmatrix_ndims;
	if (psize )*psize =outmatrix_size;
	else xDelete<int>(outmatrix_size);

}
/*}}}*/
/*FUNCTION FetchData(int** pmatrix,int* pM,int *pN,const mxArray* dataref){{{*/
void FetchData(int** pmatrix,int* pM,int *pN,const mxArray* dataref){

	int     i,outmatrix_rows,outmatrix_cols;
	double *doublematrix=NULL;
	int    *outmatrix=NULL;

	if(mxIsEmpty(dataref) ){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix_cols=0;
		outmatrix=NULL;
	}
	else if( mxIsClass(dataref,"double") || 
				mxIsClass(dataref,"single") || 
				mxIsClass(dataref,"int16") || 
				mxIsClass(dataref,"int8") || 
				mxIsClass(dataref,"uint8")){

		/*Check dataref is not pointing to NaN: */
		if ( mxIsNaN(*(mxGetPr(dataref))) && (mxGetM(dataref)==1) && (mxGetN(dataref)==1) ){
			outmatrix_rows=0;
			outmatrix_cols=0;
			outmatrix=NULL;
		}
		else{
			if(!mxIsClass(dataref,"double") && !mxIsClass(dataref,"single")){
				_printf_("Warning: converting matlab data from '" << mxGetClassName(dataref) << "' to 'double'\n");
			}
			/*Convert matlab matrix to double* matrix: */
			MatlabMatrixToDoubleMatrix(&doublematrix,&outmatrix_rows,&outmatrix_cols,dataref);

			/*Convert double matrix into integer matrix: */
			outmatrix=xNew<int>(outmatrix_rows*outmatrix_cols);
			for(i=0;i<outmatrix_rows*outmatrix_cols;i++)outmatrix[i]=(int)doublematrix[i];
		}
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM)*pM=outmatrix_rows;
	if (pN)*pN=outmatrix_cols;
}
/*}}}*/
/*FUNCTION FetchData(bool** pmatrix,int* pM,int *pN,const mxArray* dataref){{{*/
void FetchData(bool** pmatrix,int* pM,int *pN,const mxArray* dataref){

	int     i,outmatrix_rows,outmatrix_cols;
	double *doublematrix=NULL;
	bool   *outmatrix=NULL;

	if(mxIsEmpty(dataref) ){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix_cols=0;
		outmatrix=NULL;
	}
	else if (mxIsClass(dataref,"double") ){

		/*Check dataref is not pointing to NaN: */
		if ( mxIsNaN(*(mxGetPr(dataref))) && (mxGetM(dataref)==1) && (mxGetN(dataref)==1) ){
			outmatrix_rows=0;
			outmatrix_cols=0;
			outmatrix=NULL;
		}
		else{

			/*Convert matlab matrix to double* matrix: */
			MatlabMatrixToDoubleMatrix(&doublematrix,&outmatrix_rows,&outmatrix_cols,dataref);

			/*Convert double matrix into integer matrix: */
			outmatrix=xNew<bool>(outmatrix_rows*outmatrix_cols);
			for(i=0;i<outmatrix_rows;i++)outmatrix[i]=(bool)doublematrix[i];
		}
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM)*pM=outmatrix_rows;
	if (pN)*pN=outmatrix_cols;
}
/*}}}*/
/*FUNCTION FetchData(bool** pmatrix,int* pnumel,int* pndims,int** psize,const mxArray* dataref){{{*/
void FetchData(bool** pmatrix,int* pnumel,int* pndims,int** psize,const mxArray* dataref){

	int      i;
	int      outmatrix_numel,outmatrix_ndims;
	int*     outmatrix_size=NULL;
	double*  doublematrix=NULL;
	bool*    outmatrix=NULL;

	if(mxIsEmpty(dataref) ){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_numel=0;
		outmatrix_ndims=0;
		outmatrix_size =NULL;
		outmatrix=NULL;
	}
	else if (mxIsClass(dataref,"logical") ){

		/*Check dataref is not pointing to NaN: */
		if ( mxIsNaN(*((bool*)mxGetData(dataref))) && (mxGetNumberOfElements(dataref)==1) ){
			outmatrix_numel=0;
			outmatrix_ndims=0;
			outmatrix_size =NULL;
			outmatrix=NULL;
		}
		else{

			/*Convert matlab n-dim array to bool* matrix: */
			MatlabNArrayToNArray(&outmatrix,&outmatrix_numel,&outmatrix_ndims,&outmatrix_size,dataref);
		}
	}
	else if (mxIsClass(dataref,"double") ){

		/*Check dataref is not pointing to NaN: */
		if ( mxIsNaN(*(mxGetPr(dataref))) && (mxGetNumberOfElements(dataref)==1) ){
			outmatrix_numel=0;
			outmatrix_ndims=0;
			outmatrix_size =NULL;
			outmatrix=NULL;
		}
		else{

			/*Convert matlab n-dim array to double* matrix: */
			MatlabNArrayToNArray(&doublematrix,&outmatrix_numel,&outmatrix_ndims,&outmatrix_size,dataref);

			/*Convert double matrix into bool matrix: */
			outmatrix=xNew<bool>(outmatrix_numel);
			for(i=0;i<outmatrix_numel;i++)outmatrix[i]=(bool)doublematrix[i];
			xDelete<double>(doublematrix);
		}
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pnumel)*pnumel=outmatrix_numel;
	if (pndims)*pndims=outmatrix_ndims;
	if (psize )*psize =outmatrix_size;
	else xDelete<int>(outmatrix_size);

}
/*}}}*/
/*FUNCTION FetchData(double** pvector,int* pM,const mxArray* dataref){{{*/
void FetchData(double** pvector,int* pM,const mxArray* dataref){

	double* outvector=NULL;
	int M,N;

	/*Use Fetch matrix*/
	FetchData(&outvector,&M,&N,dataref) ;

	/*Check that it is a vector*/
	if(M*N>0 && (M!=1 && N!=1)){
		_error_("input vector of size " << M << "x" << N << " should have only one column");
	}

	/*Transpose Row vectors*/
	if(M==1 && N>1) M=N;

	/*Assign output pointers:*/
	*pvector=outvector;
	if(pM)*pM=M;
}
/*}}}*/
/*FUNCTION FetchData(int** pvector,int* pM,const mxArray* dataref){{{*/
void FetchData(int** pvector,int* pM,const mxArray* dataref){

	int    i;
	double *doublevector   = NULL;
	int    *outvector      = NULL;
	int     outvector_rows;

	if(mxIsEmpty(dataref)){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outvector_rows=0;
		outvector=NULL;
	}
	else if (mxIsClass(dataref,"double") ){

		/*Convert matlab vector to double*  vector: */
		FetchData(&doublevector,&outvector_rows,dataref);

		/*Convert double vector into integer vector: */
		outvector=xNew<int>(outvector_rows);
		for(i=0;i<outvector_rows;i++)outvector[i]=(int)doublevector[i];
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pvector=outvector;
	if (pM)*pM=outvector_rows;
}
/*}}}*/
/*FUNCTION FetchData(bool** pvector,int* pM,const mxArray* dataref){{{*/
void FetchData(bool** pvector,int* pM,const mxArray* dataref){

	int    i;
	double *doublevector   = NULL;
	bool   *outvector      = NULL;
	int     outvector_rows;

	if(mxIsEmpty(dataref)){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outvector_rows=0;
		outvector=NULL;
	}
	else if (mxIsClass(dataref,"double") ){

		/*Convert matlab vector to double*  vector: */
		FetchData(&doublevector,&outvector_rows,dataref);

		/*Convert double vector into integer vector: */
		outvector=xNew<bool>(outvector_rows);
		for(i=0;i<outvector_rows;i++)outvector[i]=(bool)doublevector[i];
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pvector=outvector;
	if (pM)*pM=outvector_rows;
}
/*}}}*/
/*FUNCTION FetchData(float** pvector,int* pM,const mxArray* dataref){{{*/
void FetchData(float** pvector,int* pM,const mxArray* dataref){

	int    i;
	double *doublevector   = NULL;
	float  *outvector      = NULL;
	int     outvector_rows;

	if(mxIsEmpty(dataref)){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outvector_rows=0;
		outvector=NULL;
	}
	else if (mxIsClass(dataref,"double") ){

		/*Convert matlab vector to double*  vector: */
		FetchData(&doublevector,&outvector_rows,dataref);

		/*Convert double vector into float vector: */
		outvector=xNew<float>(outvector_rows);
		for(i=0;i<outvector_rows;i++)outvector[i]=(float)doublevector[i];
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pvector=outvector;
	if (pM)*pM=outvector_rows;
}
/*}}}*/
/*FUNCTION FetchData(char** pstring,const mxArray* dataref){{{*/
void FetchData(char** pstring,const mxArray* dataref){

	char* outstring=NULL;

	/*Ok, the string should be coming directly from the matlab workspace: */
	if (!mxIsClass(dataref,"char")){
		_error_("input data_type is not a string!");
	}
	else{
		/*Recover the string:*/
		int stringlen;

		stringlen = mxGetM(dataref)*mxGetN(dataref)+1;
		outstring =xNew<char>(stringlen);
		mxGetString(dataref,outstring,stringlen);
	}

	/*Assign output pointers:*/
	*pstring=outstring;
}/*}}}*/
/*FUNCTION FetchData(char** pmatrix,int* pnumel,int* pndims,int** psize,const mxArray* dataref){{{*/
void FetchData(char** pmatrix,int* pnumel,int* pndims,int** psize,const mxArray* dataref){

	int      outmatrix_numel,outmatrix_ndims;
	int*     outmatrix_size=NULL;
	char*    outmatrix=NULL;

	if(mxIsEmpty(dataref) ){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_numel=0;
		outmatrix_ndims=0;
		outmatrix_size =NULL;
		outmatrix=NULL;
	}
	else if (mxIsClass(dataref,"char") ){

		/*Convert matlab n-dim array to char* matrix: */
		MatlabNArrayToNArray(&outmatrix,&outmatrix_numel,&outmatrix_ndims,&outmatrix_size,dataref);
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pnumel)*pnumel=outmatrix_numel;
	if (pndims)*pndims=outmatrix_ndims;
	if (psize )*psize =outmatrix_size;
	else xDelete<int>(outmatrix_size);

}
/*}}}*/
/*FUNCTION FetchData(double* pscalar,const mxArray* dataref){{{*/
void FetchData(double* pscalar,const mxArray* dataref){

	double scalar;

	if (!mxIsClass(dataref,"double")){
		_error_("input data_type is not a double!");
	}
	else{
		/*Recover the double: */
		scalar=mxGetScalar(dataref);
	}

	/*Assign output pointers:*/
	*pscalar=scalar;
}
/*}}}*/
/*FUNCTION FetchData(int* pinteger,const mxArray* dataref){{{*/
void FetchData(int* pinteger,const mxArray* dataref){

	int integer;

	if (!mxIsClass(dataref,"double")){
		_error_("input data_type is not a scalar!");
	}
	else{
		/*Recover the double: */
		integer=(int)mxGetScalar(dataref);
	}

	/*Assign output pointers:*/
	*pinteger=integer;
}
/*}}}*/
/*FUNCTION FetchData(bool* pboolean,const mxArray* dataref){{{*/
void FetchData(bool* pboolean,const mxArray* dataref){

	bool* mxbool_ptr=NULL;

	if (mxIsClass(dataref,"logical")){
		if(mxGetM(dataref)!=1) _error_("input data is not of size 1x1");
		if(mxGetN(dataref)!=1) _error_("input data is not of size 1x1");
		mxbool_ptr=mxGetLogicals(dataref);
	}
	else{
		_error_("input data_type is not a bool!");
	}

	*pboolean=*mxbool_ptr;
}
/*}}}*/

/*ISSM objects*/
/*FUNCTION FetchData(BamgGeom** pbamggeom,const mxArray* dataref){{{*/
void FetchData(BamgGeom** pbamggeom,const mxArray* dataref){

	/*Initialize output*/
	BamgGeom* bamggeom=new BamgGeom();

	/*Fetch all fields*/
	FetchData(&bamggeom->Vertices,&bamggeom->VerticesSize[0],&bamggeom->VerticesSize[1],mxGetAssignedField(dataref,0,"Vertices"));
	FetchData(&bamggeom->Edges, &bamggeom->EdgesSize[0], &bamggeom->EdgesSize[1], mxGetAssignedField(dataref,0,"Edges"));
	FetchData(&bamggeom->Corners, &bamggeom->CornersSize[0], &bamggeom->CornersSize[1], mxGetAssignedField(dataref,0,"Corners"));
	FetchData(&bamggeom->RequiredVertices,&bamggeom->RequiredVerticesSize[0],&bamggeom->RequiredVerticesSize[1],mxGetAssignedField(dataref,0,"RequiredVertices"));
	FetchData(&bamggeom->RequiredEdges, &bamggeom->RequiredEdgesSize[0], &bamggeom->RequiredEdgesSize[1], mxGetAssignedField(dataref,0,"RequiredEdges"));
	FetchData(&bamggeom->CrackedEdges,&bamggeom->CrackedEdgesSize[0],&bamggeom->CrackedEdgesSize[1],mxGetAssignedField(dataref,0,"CrackedEdges"));
	FetchData(&bamggeom->SubDomains,&bamggeom->SubDomainsSize[0],&bamggeom->SubDomainsSize[1],mxGetAssignedField(dataref,0,"SubDomains"));

	/*Assign output pointers:*/
	*pbamggeom=bamggeom;
}
/*}}}*/
/*FUNCTION FetchData(BamgMesh** pbamgmesh,const mxArray* dataref){{{*/
void FetchData(BamgMesh** pbamgmesh,const mxArray* dataref){

	/*Initialize output*/
	BamgMesh* bamgmesh=new BamgMesh();

	/*Fetch all fields*/
	FetchData(&bamgmesh->Vertices,&bamgmesh->VerticesSize[0],&bamgmesh->VerticesSize[1],mxGetAssignedField(dataref,0,"Vertices"));
	FetchData(&bamgmesh->Edges, &bamgmesh->EdgesSize[0], &bamgmesh->EdgesSize[1], mxGetAssignedField(dataref,0,"Edges"));
	FetchData(&bamgmesh->Triangles, &bamgmesh->TrianglesSize[0], &bamgmesh->TrianglesSize[1], mxGetAssignedField(dataref,0,"Triangles"));
	FetchData(&bamgmesh->CrackedEdges,&bamgmesh->CrackedEdgesSize[0],&bamgmesh->CrackedEdgesSize[1],mxGetAssignedField(dataref,0,"CrackedEdges"));
	FetchData(&bamgmesh->VerticesOnGeomEdge,&bamgmesh->VerticesOnGeomEdgeSize[0],&bamgmesh->VerticesOnGeomEdgeSize[1],mxGetAssignedField(dataref,0,"VerticesOnGeomEdge"));
	FetchData(&bamgmesh->VerticesOnGeomVertex,&bamgmesh->VerticesOnGeomVertexSize[0],&bamgmesh->VerticesOnGeomVertexSize[1],mxGetAssignedField(dataref,0,"VerticesOnGeomVertex"));
	FetchData(&bamgmesh->EdgesOnGeomEdge, &bamgmesh->EdgesOnGeomEdgeSize[0], &bamgmesh->EdgesOnGeomEdgeSize[1], mxGetAssignedField(dataref,0,"EdgesOnGeomEdge"));
	FetchData(&bamgmesh->IssmSegments,&bamgmesh->IssmSegmentsSize[0],&bamgmesh->IssmSegmentsSize[1],mxGetAssignedField(dataref,0,"IssmSegments"));

	/*Assign output pointers:*/
	*pbamgmesh=bamgmesh;
}
/*}}}*/
/*FUNCTION FetchData(BamgOpts** pbamgopts,const mxArray* dataref){{{*/
void FetchData(BamgOpts** pbamgopts,const mxArray* dataref){

	/*Initialize output*/
	BamgOpts* bamgopts=new BamgOpts();

	/*Fetch all fields*/
	FetchData(&bamgopts->anisomax,mxGetField(dataref,0,"anisomax"));
	FetchData(&bamgopts->cutoff,mxGetField(dataref,0,"cutoff"));
	FetchData(&bamgopts->coeff,mxGetField(dataref,0,"coeff"));
	FetchData(&bamgopts->errg,mxGetField(dataref,0,"errg"));
	FetchData(&bamgopts->gradation,mxGetField(dataref,0,"gradation"));
	FetchData(&bamgopts->Hessiantype,mxGetField(dataref,0,"Hessiantype"));
	FetchData(&bamgopts->MaxCornerAngle,mxGetField(dataref,0,"MaxCornerAngle"));
	FetchData(&bamgopts->maxnbv,mxGetField(dataref,0,"maxnbv"));
	FetchData(&bamgopts->maxsubdiv,mxGetField(dataref,0,"maxsubdiv"));
	FetchData(&bamgopts->Metrictype,mxGetField(dataref,0,"Metrictype"));
	FetchData(&bamgopts->nbjacobi,mxGetField(dataref,0,"nbjacobi"));
	FetchData(&bamgopts->nbsmooth,mxGetField(dataref,0,"nbsmooth"));
	FetchData(&bamgopts->omega,mxGetField(dataref,0,"omega"));
	FetchData(&bamgopts->power,mxGetField(dataref,0,"power"));
	FetchData(&bamgopts->random,mxGetField(dataref,0,"random"));
	FetchData(&bamgopts->verbose,mxGetField(dataref,0,"verbose"));

	FetchData(&bamgopts->Crack,mxGetField(dataref,0,"Crack"));
	FetchData(&bamgopts->geometricalmetric,mxGetField(dataref,0,"geometricalmetric"));
	FetchData(&bamgopts->KeepVertices,mxGetField(dataref,0,"KeepVertices"));
	FetchData(&bamgopts->splitcorners,mxGetField(dataref,0,"splitcorners"));

	FetchData(&bamgopts->hmin,mxGetField(dataref,0,"hmin"));
	FetchData(&bamgopts->hmax,mxGetField(dataref,0,"hmax"));
	FetchData(&bamgopts->hminVertices,&bamgopts->hminVerticesSize[0],&bamgopts->hminVerticesSize[1],mxGetField(dataref,0,"hminVertices"));
	FetchData(&bamgopts->hmaxVertices,&bamgopts->hmaxVerticesSize[0],&bamgopts->hmaxVerticesSize[1],mxGetField(dataref,0,"hmaxVertices"));
	FetchData(&bamgopts->hVertices,&bamgopts->hVerticesSize[0],&bamgopts->hVerticesSize[1],mxGetField(dataref,0,"hVertices"));
	FetchData(&bamgopts->metric,&bamgopts->metricSize[0],&bamgopts->metricSize[1],mxGetField(dataref,0,"metric"));
	FetchData(&bamgopts->field,&bamgopts->fieldSize[0],&bamgopts->fieldSize[1],mxGetField(dataref,0,"field"));
	FetchData(&bamgopts->err,&bamgopts->errSize[0],&bamgopts->errSize[1],mxGetField(dataref,0,"err"));

	/*Additional checks*/
	bamgopts->Check();

	/*Assign output pointers:*/
	*pbamgopts=bamgopts;
}
/*}}}*/
/*FUNCTION FetchData(Options** poptions,const mxArray** pdataref){{{*/
void FetchData(Options** poptions,int istart, int nrhs,const mxArray** pdataref){

	char   *name   = NULL;
	Option *option = NULL;

	/*Initialize output*/
	Options* options=new Options();

	/*Fetch all options*/
	for (int i=istart; i<nrhs; i=i+2){
		if (!mxIsClass(pdataref[i],"char")) _error_("Argument " << i+1 << " must be name of option");

		FetchData(&name,pdataref[i]);
		if(i+1 == nrhs) _error_("Argument " << i+2 << " must exist and be value of option \"" << name << "\".");

		option=(Option*)OptionParse(name,&pdataref[i+1]);
		options->AddOption(option);
		option=NULL;
	}

	/*Assign output pointers:*/
	*poptions=options;
}
/*}}}*/
/*FUNCTION FetchData(Contours** pcontours,const mxArray* dataref){{{*/
void FetchData(Contours** pcontours,const mxArray* dataref){

	int             numcontours,index,test1,test2;
	char            *contourname = NULL;
	Contours        *contours    = NULL;
	Contour<double> *contouri    = NULL;

	if (mxIsClass(dataref,"char")){
		FetchData(&contourname,dataref);
		contours=ExpRead<double>(contourname);
	}
	else if(mxIsClass(dataref,"struct")){

		contours=new Contours();
		numcontours=mxGetNumberOfElements(dataref);

		for(int i=0;i<numcontours;i++){

			contouri=new Contour<double>();

			index = mxGetFieldNumber(dataref,"nods");
			if(index==-1) _error_("input structure does not have a 'nods' field");
			FetchData(&contouri->nods,mxGetFieldByNumber(dataref,i,index));

			index = mxGetFieldNumber(dataref,"x");
			if(index==-1) _error_("input structure does not have a 'x' field");
			FetchData(&contouri->x,&test1,&test2,mxGetFieldByNumber(dataref,i,index));
			if(test1!=contouri->nods || test2!=1) _error_("field x should be of size ["<<contouri->nods<<" 1]");

			index = mxGetFieldNumber(dataref,"y");
			if(index==-1) _error_("input structure does not have a 'y' field");
			FetchData(&contouri->y,&test1,&test2,mxGetFieldByNumber(dataref,i,index));
			if(test1!=contouri->nods || test2!=1) _error_("field y should be of size ["<<contouri->nods<<" 1]");

			contours->AddObject(contouri);
		}
	}
	else{
		_error_("Contour is neither a string nor a structure and cannot be loaded ("<<mxGetClassName(dataref)<<" not supported)");
	}

	/*clean-up and assign output pointer*/
	xDelete<char>(contourname);
	*pcontours=contours;
}
/*}}}*/

/*Toolkit*/
/*FUNCTION MatlabMatrixToDoubleMatrix {{{*/
int MatlabMatrixToDoubleMatrix(double** pmatrix,int* pmatrix_rows,int* pmatrix_cols,const mxArray* mxmatrix){

	int        i,j,count,rows,cols;

	/*output: */
	double* matrix=NULL;

	/*matlab indices: */
	mwIndex*    ir=NULL;
	mwIndex*    jc=NULL;

	/*Ok, first check if we are dealing with a sparse or full matrix: */
	if (mxIsSparse(mxmatrix)){

		/*Dealing with sparse matrix: recover size first: */
		double* pmxmatrix=(double*)mxGetPr(mxmatrix);
		rows=mxGetM(mxmatrix);
		cols=mxGetN(mxmatrix);

		if(rows*cols){
			matrix=xNewZeroInit<double>(rows*cols);

			/*Now, get ir,jc and pr: */
			ir=mxGetIr(mxmatrix);
			jc=mxGetJc(mxmatrix);

			/*Now, start inserting data into double* matrix: */
			count=0;
			for(i=0;i<cols;i++){
				for(j=0;j<(jc[i+1]-jc[i]);j++){
					matrix[rows*ir[count]+i]=pmxmatrix[count];
					count++;
				}
			}
		}

	}
	else if(mxIsClass(mxmatrix,"double")){
		/*Dealing with dense matrix: recover pointer and size: */
		double* pmxmatrix=(double*)mxGetPr(mxmatrix);
		rows=mxGetM(mxmatrix);
		cols=mxGetN(mxmatrix);

		/*Create serial matrix: */
		if(rows*cols){
			matrix=xNewZeroInit<double>(rows*cols);

			for(i=0;i<rows;i++){
				for(j=0;j<cols;j++){
					matrix[cols*i+j]=(double)pmxmatrix[rows*j+i];
				}
			}
		}
	}
	else if(mxIsClass(mxmatrix,"single")){
		/*Dealing with dense matrix: recover pointer and size: */
		float *pmxmatrix=(float*)mxGetPr(mxmatrix);
		rows=mxGetM(mxmatrix);
		cols=mxGetN(mxmatrix);

		/*Create serial matrix: */
		if(rows*cols){
			matrix=xNewZeroInit<double>(rows*cols);

			for(i=0;i<rows;i++){
				for(j=0;j<cols;j++){
					matrix[cols*i+j]=(double)pmxmatrix[rows*j+i];
				}
			}
		}
	}
	else if(mxIsClass(mxmatrix,"int16")){
		/*Dealing with dense matrix: recover pointer and size: */
		short int *pmxmatrix=(short*)mxGetPr(mxmatrix);
		rows=mxGetM(mxmatrix);
		cols=mxGetN(mxmatrix);

		/*Create serial matrix: */
		if(rows*cols){
			matrix=xNewZeroInit<double>(rows*cols);

			for(i=0;i<rows;i++){
				for(j=0;j<cols;j++){
					matrix[cols*i+j]=(double)pmxmatrix[rows*j+i];
				}
			}
		}
	}
	else if(mxIsClass(mxmatrix,"uint8")){
		/*Dealing with dense matrix: recover pointer and size: */
		char *pmxmatrix=(char*)mxGetPr(mxmatrix);
		rows=mxGetM(mxmatrix);
		cols=mxGetN(mxmatrix);

		/*Create serial matrix: */
		if(rows*cols){
			matrix=xNewZeroInit<double>(rows*cols);

			for(i=0;i<rows;i++){
				for(j=0;j<cols;j++){
					matrix[cols*i+j]=(double)pmxmatrix[rows*j+i];
				}
			}
		}
	}
	else{
		_error_("Matlab matrix type Not implemented yet");
	}

	/*Assign output pointer: */
	*pmatrix=matrix;
	*pmatrix_rows=rows;
	*pmatrix_cols=cols;

	return 1;
}/*}}}*/
/*FUNCTION MatlabNArrayToNArray(double** pmatrix,int* pmatrix_numel,int* pmatrix_ndims,int** pmatrix_size,const mxArray* mxmatrix){{{*/
int MatlabNArrayToNArray(double** pmatrix,int* pmatrix_numel,int* pmatrix_ndims,int** pmatrix_size,const mxArray* mxmatrix){

	int  i,j,rows,cols;
	int  numel,ndims;
	int *size,*dims;
	double* mxmatrix_ptr=NULL;
	const mwSize* ipt=NULL;

	/*output: */
	double* matrix=NULL;

	/*matlab indices: */
	mwIndex *ir    = NULL;
	mwIndex *jc    = NULL;
	double  *pr    = NULL;
	int      count;

	/*get Matlab matrix information: */
	numel=mxGetNumberOfElements(mxmatrix);
	ndims=mxGetNumberOfDimensions(mxmatrix);
	ipt  =mxGetDimensions(mxmatrix);
	size =xNew<int>(ndims);
	for (i=0;i<ndims;i++) size[i]=(int)ipt[i];

	/*Ok, first check if we are dealing with a sparse or full matrix: */
	if (mxIsSparse(mxmatrix)){

		/*Dealing with sparse matrix: recover size first: */
		rows = mxGetM(mxmatrix);
		cols = mxGetN(mxmatrix);

		matrix=xNewZeroInit<double>(rows*cols);

		/*Now, get ir,jc and pr: */
		ir = mxGetIr(mxmatrix);
		jc = mxGetJc(mxmatrix);
		pr = mxGetPr(mxmatrix);

		/*Now, start inserting data into double* matrix: */
		count=0;
		for(i=0;i<cols;i++){
			for(j=0;j<(jc[i+1]-jc[i]);j++){
				*(matrix+rows*ir[count]+i)=pr[count];
				count++;
			}
		}

	}
	else{

		/*Dealing with dense matrix: recover pointer and size: */
		mxmatrix_ptr=(double*)mxGetPr(mxmatrix);

		/*Create serial matrix: */
		matrix=xNewZeroInit<double>(numel);

		dims=xNew<int>(ndims);
		for(i=0;i<numel;i++){
			ColumnWiseDimsFromIndex(dims,i,size,ndims);
			j = IndexFromRowWiseDims(dims,size,ndims);
			matrix[j]=(double)mxmatrix_ptr[i];
		}
		xDelete<int>(dims);
	}

	/*Assign output pointer: */
	*pmatrix       = matrix;
	*pmatrix_numel = numel;
	*pmatrix_ndims = ndims;
	*pmatrix_size  = size;

	return 1;
}
/*}}}*/
/*FUNCTION MatlabNArrayToNArray(bool** pmatrix,int* pmatrix_numel,int* pmatrix_ndims,int** pmatrix_size,const mxArray* mxmatrix){{{*/
int MatlabNArrayToNArray(bool** pmatrix,int* pmatrix_numel,int* pmatrix_ndims,int** pmatrix_size,const mxArray* mxmatrix){

	int  i,j,rows,cols;
	int  numel,ndims;
	int *size,*dims;
	bool* mxmatrix_ptr=NULL;
	const mwSize* ipt=NULL;

	/*output: */
	bool* matrix=NULL;

	/*matlab indices: */
	mwIndex *ir    = NULL;
	mwIndex *jc    = NULL;
	bool    *pm    = NULL;
	int      count;

	/*get Matlab matrix information: */
	numel = mxGetNumberOfElements(mxmatrix);
	ndims = mxGetNumberOfDimensions(mxmatrix);
	ipt   = mxGetDimensions(mxmatrix);
	size  = xNew<int>(ndims);
	for (i=0;i<ndims;i++) size[i]=(int)ipt[i];

	/*Ok, first check if we are dealing with a sparse or full matrix: */
	if (mxIsSparse(mxmatrix)){

		/*Dealing with sparse matrix: recover size first: */
		rows=mxGetM(mxmatrix);
		cols=mxGetN(mxmatrix);
		matrix=xNewZeroInit<bool>(rows*cols);

		/*Now, get ir,jc and pm: */
		ir=mxGetIr(mxmatrix);
		jc=mxGetJc(mxmatrix);
		pm=(bool*)mxGetData(mxmatrix);

		/*Now, start inserting data into bool* matrix: */
		count=0;
		for(i=0;i<cols;i++){
			for(j=0;j<(jc[i+1]-jc[i]);j++){
				matrix[rows*ir[count]+i]=pm[count];
				count++;
			}
		}
	}
	else{

		/*Dealing with dense matrix: recover pointer and size: */
		mxmatrix_ptr=(bool*)mxGetData(mxmatrix);

		/*Create serial matrix: */
		matrix=xNew<bool>(numel);
		dims=xNew<int>(ndims);
		for(i=0;i<numel;i++){
			ColumnWiseDimsFromIndex(dims,i,size,ndims);
			j=IndexFromRowWiseDims(dims,size,ndims);
			matrix[j]=(bool)mxmatrix_ptr[i];
		}
		xDelete<int>(dims);
	}

	/*Assign output pointer: */
	*pmatrix       = matrix;
	*pmatrix_numel = numel;
	*pmatrix_ndims = ndims;
	*pmatrix_size  = size;

	return 1;
}
/*}}}*/
/*FUNCTION MatlabNArrayToNArray(char** pmatrix,int* pmatrix_numel,int* pmatrix_ndims,int** pmatrix_size,const mxArray* mxmatrix){{{*/
int MatlabNArrayToNArray(char** pmatrix,int* pmatrix_numel,int* pmatrix_ndims,int** pmatrix_size,const mxArray* mxmatrix){

	int           i,j,rows,cols;
	int           numel,ndims;
	int          *size , *dims;
	mxChar       *mxmatrix_ptr = NULL;
	const mwSize *ipt          = NULL;

	/*output: */
	char* matrix=NULL;

	/*matlab indices: */
	mwIndex *ir    = NULL;
	mwIndex *jc    = NULL;
	char    *pm    = NULL;
	int      count;

	/*get Matlab matrix information: */
	numel = mxGetNumberOfElements(mxmatrix);
	ndims = mxGetNumberOfDimensions(mxmatrix);
	ipt   = mxGetDimensions(mxmatrix);
	size  = xNew<int>(ndims);
	for (i=0;i<ndims;i++) size[i]=(int)ipt[i];

	/*Ok, first check if we are dealing with a sparse or full matrix: */
	if (mxIsSparse(mxmatrix)){

		/*Dealing with sparse matrix: recover size first: */
		rows = mxGetM(mxmatrix);
		cols = mxGetN(mxmatrix);
		matrix=xNew<char>(rows*cols);

		/*Now, get ir,jc and pm: */
		ir = mxGetIr(mxmatrix);
		jc = mxGetJc(mxmatrix);
		pm = (char*)mxGetData(mxmatrix);

		/*Now, start inserting data into char* matrix: */
		count=0;
		for(i=0;i<cols;i++){
			for(j=0;j<(jc[i+1]-jc[i]);j++){
				matrix[rows*ir[count]+i]=(char)pm[count];
				count++;
			}
		}
	}
	else{
		/*Dealing with dense matrix: recover pointer and size: */
		mxmatrix_ptr=mxGetChars(mxmatrix);

		/*Create serial matrix: */
		matrix=xNew<char>(numel+1);
		matrix[numel]='\0';

		/*looping code adapted from Matlab example explore.c: */
		int elements_per_page = size[0] * size[1];
		/* total_number_of_pages = size[2] x size[3] x ... x size[N-1] */
		int total_number_of_pages = 1;
		for (i=2; i<ndims; i++) {
			total_number_of_pages *= size[i];
		}

		i=0;
		for (int page=0; page < total_number_of_pages; page++) {
			int row;
			/* On each page, walk through each row. */
			for (row=0; row<size[0]; row++)  {
				int column;
				j = (page * elements_per_page) + row;

				/* Walk along each column in the current row. */
				for (column=0; column<size[1]; column++) {
					*(matrix+i++)=(char)*(mxmatrix_ptr+j);
					j += size[0];
				}
			}
		}
	}

	/*Assign output pointer: */
	*pmatrix       = matrix;
	*pmatrix_numel = numel;
	*pmatrix_ndims = ndims;
	*pmatrix_size  = size;

	return 1;
}
/*}}}*/
/*FUNCTION mxGetAssignedField{{{*/
mxArray* mxGetAssignedField(const mxArray* pmxa_array,int number,const char* field){

	//output
	mxArray* mxfield=NULL;

	//input
	mxArray    *inputs[2];
	mxArray    *pindex      = NULL;
	const char *fnames[2];
	mwSize      ndim        = 2;
	mwSize      onebyone[2] = {1,1};

	//We want to call the subsasgn method, and get the returned array.This ensures that if we are running 
	//large sized problems, the data is truly loaded from disk by the model subsasgn class method.
	inputs[0]=(mxArray*)pmxa_array; //this is the model

	//create index structure used in the assignment (index.type='.' and index.subs='x' for field x for ex)
	fnames[0] = "type";
	fnames[1] = "subs";
	pindex=mxCreateStructArray( ndim,onebyone,2,fnames);
	mxSetField( pindex, 0, "type",mxCreateString("."));
	mxSetField( pindex, 0, "subs",mxCreateString(field));
	inputs[1]=pindex;

	mexCallMATLAB( 1, &mxfield, 2, (mxArray**)inputs, "subsref");

	return mxfield;
}/*}}}*/

GenericOption<double>* OptionDoubleParse( char* name, const mxArray* prhs[]){ /*{{{*/

	GenericOption<double> *odouble = NULL;

	/*check and parse the name  */
	odouble=new GenericOption<double>();
	odouble->name =xNew<char>(strlen(name)+1);
	memcpy(odouble->name,name,(strlen(name)+1)*sizeof(char));
	FetchData(&odouble->value,prhs[0]);
	odouble->numel=1;
	odouble->ndims=1;
	odouble->size=NULL;

	return(odouble);
}/*}}}*/
GenericOption<double*>* OptionDoubleArrayParse( char* name, const mxArray* prhs[]){ /*{{{*/

	GenericOption<double*> *odouble = NULL;

	/*check and parse the name  */
	odouble=new GenericOption<double*>();
	odouble->name =xNew<char>(strlen(name)+1);
	memcpy(odouble->name,name,(strlen(name)+1)*sizeof(char));

	/*check and parse the value  */
	if (!mxIsClass(prhs[0],"double")){
		_error_("Value of option \"" << odouble->name  << "\" must be class \"double\", not class \"" << mxGetClassName(prhs[0]) <<"\".");
	}
	FetchData(&odouble->value,&odouble->numel,&odouble->ndims,&odouble->size,prhs[0]);

	return(odouble);
}/*}}}*/
GenericOption<bool*>* OptionLogicalParse( char* name, const mxArray* prhs[]){ /*{{{*/

	GenericOption<bool*> *ological = NULL;

	/*check and parse the name  */
	ological=new GenericOption<bool*>();
	ological->name =xNew<char>(strlen(name)+1);
	memcpy(ological->name,name,(strlen(name)+1)*sizeof(char));

	/*check and parse the value  */
	if (!mxIsClass(prhs[0],"logical")){
		_error_("Value of option \"" << ological->name  << "\" must be class \"logical\", not class \"" << mxGetClassName(prhs[0]) <<"\".");
	}
	FetchData(&ological->value,&ological->numel,&ological->ndims,&ological->size,prhs[0]);

	return(ological);
}/*}}}*/
GenericOption<char*>* OptionCharParse( char* name, const mxArray* prhs[]){ /*{{{*/

	GenericOption<char*>  *ochar = NULL;

	/*check and parse the name  */
	ochar=new GenericOption<char*>();
	ochar->name =xNew<char>(strlen(name)+1);
	memcpy(ochar->name,name,(strlen(name)+1)*sizeof(char));

	/*check and parse the value  */
	if (!mxIsClass(prhs[0],"char")){
		_error_("Value of option \"" << ochar->name  << "\" must be class \"char\", not class \"" << mxGetClassName(prhs[0]) <<"\".");
	}
	FetchData(&ochar->value,&ochar->numel,&ochar->ndims,&ochar->size,prhs[0]);

	return(ochar);
}/*}}}*/
GenericOption<Options**>* OptionStructParse( char* name, const mxArray* prhs[]){ /*{{{*/

	int            i;
	char           namei[161];
	Option*                   option      = NULL;
	GenericOption<Options**>  *ostruct    = NULL;
	const mwSize  *ipt        = NULL;
	const mxArray *structi;
	mwIndex        sindex;

	/*check and parse the name  */
	ostruct=new GenericOption<Options**>();
	ostruct->name =xNew<char>(strlen(name)+1);
	memcpy(ostruct->name,name,(strlen(name)+1)*sizeof(char));

	/*check and parse the value  */
	if (!mxIsClass(prhs[0],"struct")){
		_error_("Value of option \"" << ostruct->name  << "\" must be class \"struct\", not class \"" << mxGetClassName(prhs[0]) <<"\".");
	}
	ostruct->numel=mxGetNumberOfElements(prhs[0]);
	ostruct->ndims=mxGetNumberOfDimensions(prhs[0]);
	ipt           =mxGetDimensions(prhs[0]);
	ostruct->size =xNew<int>(ostruct->ndims);
	for (i=0; i<ostruct->ndims; i++) ostruct->size[i]=(int)ipt[i];
	if (ostruct->numel) ostruct->value=xNew<Options*>(ostruct->numel);

	/*loop through and process each element of the struct array  */
	for (sindex=0; sindex<ostruct->numel; sindex++) {
		ostruct->value[sindex]=new Options;

		/*loop through and process each field for the element  */
		for (i=0; i<mxGetNumberOfFields(prhs[0]); i++) {
			sprintf(namei,"%s.%s",name,mxGetFieldNameByNumber(prhs[0],i));
			structi=mxGetFieldByNumber(prhs[0],sindex,i);

			option=(Option*)OptionParse(namei,&structi);
			ostruct->value[sindex]->AddObject((Object*)option);
			option=NULL;
		}
	}

	return(ostruct);
}/*}}}*/
GenericOption<Options*>* OptionCellParse( char* name, const mxArray* prhs[]){ /*{{{*/

	int            i;
	int           *dims;
	char           namei[161];
	char           cstr[81];
	GenericOption<Options*> *ocell      = NULL;
	Option        *option     = NULL;
	const mwSize  *ipt        = NULL;
	const mxArray *celli;
	mwIndex        cindex;

	/*check and parse the name  */
	ocell=new GenericOption<Options*>();
	ocell->name =xNew<char>(strlen(name)+1);
	memcpy(ocell->name,name,(strlen(name)+1)*sizeof(char));

	/*check and parse the value  */
	if (!mxIsClass(prhs[0],"cell")){
		_error_("Value of option \"" << ocell->name  << "\" must be class \"cell\", not class \"" << mxGetClassName(prhs[0]) <<"\".");
	}

	ocell->numel=mxGetNumberOfElements(prhs[0]);
	ocell->ndims=mxGetNumberOfDimensions(prhs[0]);
	ipt         =mxGetDimensions(prhs[0]);
	ocell->size =xNew<int>(ocell->ndims);
	for (i=0; i<ocell->ndims; i++) ocell->size[i]=(int)ipt[i];
	ocell->value=new Options;

	/*loop through and process each element of the cell array  */
	dims=xNew<int>(ocell->ndims);
	for (cindex=0; cindex<ocell->numel; cindex++) {
		ColumnWiseDimsFromIndex(dims,(int)cindex,ocell->size,ocell->ndims);
		StringFromDims(cstr,dims,ocell->ndims);
		#ifdef _INTEL_WIN_
			_snprintf(namei,161,"%s%s",name,cstr);
		#else
			snprintf(namei,161,"%s%s",name,cstr);
		#endif
		celli=mxGetCell(prhs[0],cindex);

		option=(Option*)OptionParse(namei,&celli);
		ocell->value->AddObject((Object*)option);
		option=NULL;
	}
	xDelete<int>(dims);

	return(ocell);
}/*}}}*/
Option* OptionParse(char* name, const mxArray* prhs[]){ /*{{{*/

	Option  *option = NULL;
	mxArray *lhs[1];

	/*parse the value according to the matlab data type  */
	if     (mxIsClass(prhs[0],"double")  && (mxGetNumberOfElements(prhs[0])==1))
	 option=(Option*)OptionDoubleParse(name,prhs);
	else if(mxIsClass(prhs[0],"double")  && (mxGetNumberOfElements(prhs[0])!=1))
	 option=(Option*)OptionDoubleArrayParse(name,prhs);
	else if(mxIsClass(prhs[0],"logical"))
	 option=(Option*)OptionLogicalParse(name,prhs);
	else if(mxIsClass(prhs[0],"char"))
	 option=(Option*)OptionCharParse(name,prhs);
	else if(mxIsClass(prhs[0],"struct"))
	 option=(Option*)OptionStructParse(name,prhs);
	else if(mxIsClass(prhs[0],"cell"))
	 option=(Option*)OptionCellParse(name,prhs);
	else {
		_printf0_("  Converting value of option \"" << name << "\" from unrecognized class \"" << mxGetClassName(prhs[0]) << "\" to class \"" << "struct" << "\".\n");
		if (!mexCallMATLAB(1,lhs,1,(mxArray**)prhs,"struct")) {
			option=(Option*)OptionStructParse(name,(const mxArray**)lhs);
			mxDestroyArray(lhs[0]);
		}
		else _error_("Second argument value of option \""<< name <<"\" is of unrecognized class \""<< mxGetClassName(prhs[0]) <<"\".");
	}

	return(option);
}/*}}}*/
