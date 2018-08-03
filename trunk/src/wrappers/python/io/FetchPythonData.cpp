/*\file FetchData.cpp:
 * \brief: general I/O interface to fetch data in python
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#define PY_ARRAY_UNIQUE_SYMBOL PythonIOSymbol
#define NO_IMPORT

#include "./pythonio.h"
#include "../../c/shared/shared.h"

/*Primitive data types*/
/*FUNCTION FetchData(double* pscalar,PyObject* py_float){{{*/
void FetchData(double* pscalar,PyObject* py_float){

	double dscalar;

	/*return internal value: */
	if      (PyFloat_Check(py_float))
		dscalar=PyFloat_AsDouble(py_float);
	else if (PyLong_Check(py_float)){
		#if _PYTHON_MAJOR_ == 3
		dscalar=(double)PyLong_AsLong(py_float);
		#else
		dscalar=PyLong_AsDouble(py_float);
		#endif
	}
	else if (PyInt_Check(py_float))
	 dscalar=(double)PyInt_AsLong(py_float);
	else if (PyBool_Check(py_float))
		dscalar=(double)PyLong_AsLong(py_float);
	else if (PyTuple_Check(py_float) && (int)PyTuple_Size(py_float)==1)
		FetchData(&dscalar,PyTuple_GetItem(py_float,(Py_ssize_t)0));
	else if (PyList_Check(py_float) && (int)PyList_Size(py_float)==1)
		FetchData(&dscalar,PyList_GetItem(py_float,(Py_ssize_t)0));
	else
		_error_("unrecognized float type in input!");

	/*output: */
	*pscalar=dscalar;
}
/*}}}*/
/*FUNCTION FetchData(int* pscalar,PyObject* py_long){{{*/
void FetchData(int* pscalar, PyObject* py_long){

	int iscalar;

	/*return internal value: */
	if      (PyLong_Check(py_long))
		iscalar=(int)PyLong_AsLong(py_long);
	else if (PyInt_Check(py_long))
		iscalar=(int)PyInt_AsLong(py_long);
	else if (PyFloat_Check(py_long))
		iscalar=(int)PyFloat_AsDouble(py_long);
	else if (PyBool_Check(py_long))
		iscalar=(int)PyLong_AsLong(py_long);
	else if (PyTuple_Check(py_long) && (int)PyTuple_Size(py_long)==1)
		FetchData(&iscalar,PyTuple_GetItem(py_long,(Py_ssize_t)0));
	else if (PyList_Check(py_long) && (int)PyList_Size(py_long)==1)
		FetchData(&iscalar,PyList_GetItem(py_long,(Py_ssize_t)0));
	else
		_error_("unrecognized long type in input!");

	/*output: */
	*pscalar=iscalar;
}
/*}}}*/
/*FUNCTION FetchData(bool* pscalar,PyObject* py_boolean){{{*/
void FetchData(bool* pscalar,PyObject* py_boolean){

	bool bscalar;

	/*return internal value: */
	if      (PyBool_Check(py_boolean))
		bscalar=(bool)PyLong_AsLong(py_boolean);
	else if (PyLong_Check(py_boolean))
		bscalar=(bool)PyLong_AsLong(py_boolean);
	else if (PyLong_Check(py_boolean))
		bscalar=(bool)PyLong_AsLong(py_boolean);
	else if (PyInt_Check(py_boolean))
		bscalar=(bool)PyInt_AsLong(py_boolean);
	else if (PyTuple_Check(py_boolean) && (int)PyTuple_Size(py_boolean)==1)
		FetchData(&bscalar,PyTuple_GetItem(py_boolean,(Py_ssize_t)0));
	else if (PyList_Check(py_boolean) && (int)PyList_Size(py_boolean)==1)
		FetchData(&bscalar,PyList_GetItem(py_boolean,(Py_ssize_t)0));
	else
		_error_("unrecognized boolean type in input!");

	/*output: */
	*pscalar=bscalar;
}
/*}}}*/
/*FUNCTION FetchData(double** pmatrix,int* pM, int* pN, PyObject* py_matrix){{{*/
void FetchData(double** pmatrix,int* pM,int *pN,PyObject* py_matrix){

	/*output: */
	double* dmatrix=NULL;
	double* matrix=NULL;
	int M,N;
	int ndim;
	npy_intp*  dims=NULL;

	/*intermediary:*/
	long* lmatrix=NULL;
	bool* bmatrix=NULL;
	int i;
	PyObject* py_matrix2=NULL;

	if     (PyArray_Check((PyArrayObject*)py_matrix)) {
		/*retrieve dimensions: */
		ndim=PyArray_NDIM((const PyArrayObject*)py_matrix);
		if      (ndim==2) {
			dims=PyArray_DIMS((PyArrayObject*)py_matrix);
			M=dims[0]; N=dims[1];
		}
		else if (ndim==1) {
			dims=PyArray_DIMS((PyArrayObject*)py_matrix);
			M=dims[0]; N=1;
		}
		else
			_error_("expecting an MxN matrix or M vector in input!");

		if (M && N) {
			if (!PyArray_ISCONTIGUOUS((PyArrayObject*)py_matrix)) {
				py_matrix2=PyArray_ContiguousFromAny(py_matrix,NPY_DOUBLE,ndim,ndim);
				py_matrix=py_matrix2;
			}

			if      (PyArray_TYPE((PyArrayObject*)py_matrix) == NPY_DOUBLE) {
				/*retrieve internal value: */
				dmatrix=(double*)PyArray_DATA((PyArrayObject*)py_matrix);

				/*copy matrix: */
				matrix=xNew<double>(M*N);
//				if (PyArray_ISCONTIGUOUS((PyArrayObject*)py_matrix)) {
					memcpy(matrix,dmatrix,(M*N)*sizeof(double));
//				}

//				else {
//					int j,k,ipt=0;
//					int mstride,nstride;
//					mstride=(int)PyArray_STRIDE((PyArrayObject*)py_matrix,0)/PyArray_ITEMSIZE((PyArrayObject*)py_matrix);
//					if (ndim > 1)
//						nstride=(int)PyArray_STRIDE((PyArrayObject*)py_matrix,1)/PyArray_ITEMSIZE((PyArrayObject*)py_matrix);
//					else
//						nstride=1;
//					for (i=0; i<M; i++) {
//						k=i*mstride;
//						for (j=0; j<N; j++) {
//							matrix[ipt++]=dmatrix[k];
//							k+=nstride;
//						}
//					}
//				}
			
			}

			else if (PyArray_TYPE((PyArrayObject*)py_matrix) == NPY_LONG) {
				/*retrieve internal value: */
				lmatrix=(long*)PyArray_DATA((PyArrayObject*)py_matrix);

				/*transform into double matrix: */
				matrix=xNew<double>(M*N);
				for(i=0;i<M*N;i++)matrix[i]=(double)lmatrix[i];
			}

			else if (PyArray_TYPE((PyArrayObject*)py_matrix) == NPY_BOOL) {
				/*retrieve internal value: */
				bmatrix=(bool*)PyArray_DATA((PyArrayObject*)py_matrix);

				/*transform into double matrix: */
				matrix=xNew<double>(M*N);
				for(i=0;i<M*N;i++)matrix[i]=(double)bmatrix[i];
			}

			else
				_error_("unrecognized double pyarray type in input!");

			//if (py_matrix2) delete(py_matrix2); Not doing this for now, seems to  be creating a segfault!
		}
		else
			matrix=NULL;
	}

	else {
		M=1;
		N=1;
		matrix=xNew<double>(M*N);
		FetchData(&(matrix[0]),py_matrix);
	}

	/*output: */
	if(pM)*pM=M;
	if(pN)*pN=N;
	if(pmatrix)*pmatrix=matrix;
}
/*}}}*/
/*FUNCTION FetchData(int** pmatrix,int* pM, int* pN, PyObject* py_matrix){{{*/
void FetchData(int** pmatrix,int* pM,int *pN,PyObject* py_matrix){

	/*output: */
	int* matrix=NULL;
	int M,N;
	int ndim;
	npy_intp*  dims=NULL;

	/*intermediary:*/
	double* dmatrix=NULL;
	long* lmatrix=NULL;
	bool* bmatrix=NULL;
	int i;
	PyObject* py_matrix2=NULL;

	if     (PyArray_Check((PyArrayObject*)py_matrix)) {
		/*retrieve dimensions: */
		ndim=PyArray_NDIM((const PyArrayObject*)py_matrix);
		if      (ndim==2) {
			dims=PyArray_DIMS((PyArrayObject*)py_matrix);
			M=dims[0]; N=dims[1];
		}
		else if (ndim==1) {
			dims=PyArray_DIMS((PyArrayObject*)py_matrix);
			M=dims[0]; N=1;
		}
		else
			_error_("expecting an MxN matrix or M vector in input!");

		if (M && N) {
			if (!PyArray_ISCONTIGUOUS((PyArrayObject*)py_matrix)) {
				py_matrix2=PyArray_ContiguousFromAny(py_matrix,NPY_LONG,ndim,ndim);
				py_matrix=py_matrix2;
			}

			if      (PyArray_TYPE((PyArrayObject*)py_matrix) == NPY_DOUBLE) {
				/*retrieve internal value: */
				dmatrix=(double*)PyArray_DATA((PyArrayObject*)py_matrix);

				/*transform into integer matrix: */
				matrix=xNew<int>(M*N);
				for(i=0;i<M*N;i++)matrix[i]=(int)dmatrix[i];
			}

			else if (PyArray_TYPE((PyArrayObject*)py_matrix) == NPY_LONG) {
				/*retrieve internal value: */
				lmatrix=(long*)PyArray_DATA((PyArrayObject*)py_matrix);

				/*transform into integer matrix: */
				matrix=xNew<int>(M*N);
				for(i=0;i<M*N;i++)matrix[i]=(int)lmatrix[i];
			}

			else if (PyArray_TYPE((PyArrayObject*)py_matrix) == NPY_BOOL) {
				/*retrieve internal value: */
				bmatrix=(bool*)PyArray_DATA((PyArrayObject*)py_matrix);

				/*transform into integer matrix: */
				matrix=xNew<int>(M*N);
				for(i=0;i<M*N;i++)matrix[i]=(int)bmatrix[i];
			}

			else
				_error_("unrecognized int pyarray type in input!");

			/* These lines are causing a segfault
			if (py_matrix2)
				delete(py_matrix2);
			*/
		}
		else
			matrix=NULL;
	}

	else {
		M=1;
		N=1;
		matrix=xNew<int>(M*N);
		FetchData(&(matrix[0]),py_matrix);
	}

	/*output: */
	if(pM)*pM=M;
	if(pN)*pN=N;
	if(pmatrix)*pmatrix=matrix;
}
/*}}}*/
/*FUNCTION FetchData(bool** pmatrix,int* pM, int* pN, PyObject* py_matrix){{{*/
void FetchData(bool** pmatrix,int* pM,int *pN,PyObject* py_matrix){

	/*output: */
	bool* bmatrix=NULL;
	bool* matrix=NULL;
	int M,N;
	int ndim;
	npy_intp*  dims=NULL;

	/*intermediary:*/
	double* dmatrix=NULL;
	long* lmatrix=NULL;
	int i;
	PyObject* py_matrix2=NULL;

	if     (PyArray_Check((PyArrayObject*)py_matrix)) {
		/*retrieve dimensions: */
		ndim=PyArray_NDIM((const PyArrayObject*)py_matrix);
		if      (ndim==2) {
			dims=PyArray_DIMS((PyArrayObject*)py_matrix);
			M=dims[0]; N=dims[1];
		}
		else if (ndim==1) {
			dims=PyArray_DIMS((PyArrayObject*)py_matrix);
			M=dims[0]; N=1;
		}
		else
			_error_("expecting an MxN matrix or M vector in input!");

		if (M && N) {
			if (!PyArray_ISCONTIGUOUS((PyArrayObject*)py_matrix)) {
				py_matrix2=PyArray_ContiguousFromAny(py_matrix,NPY_BOOL,ndim,ndim);
				py_matrix=py_matrix2;
			}

			if      (PyArray_TYPE((PyArrayObject*)py_matrix) == NPY_DOUBLE) {
				/*retrieve internal value: */
				dmatrix=(double*)PyArray_DATA((PyArrayObject*)py_matrix);

				/*transform into bool matrix: */
				matrix=xNew<bool>(M*N);
				for(i=0;i<M*N;i++)matrix[i]=(bool)dmatrix[i];
			}

			else if (PyArray_TYPE((PyArrayObject*)py_matrix) == NPY_LONG) {
				/*retrieve internal value: */
				lmatrix=(long*)PyArray_DATA((PyArrayObject*)py_matrix);

				/*transform into bool matrix: */
				matrix=xNew<bool>(M*N);
				for(i=0;i<M*N;i++)matrix[i]=(bool)lmatrix[i];
			}

			else if (PyArray_TYPE((PyArrayObject*)py_matrix) == NPY_BOOL) {
				/*retrieve internal value: */
				bmatrix=(bool*)PyArray_DATA((PyArrayObject*)py_matrix);

				/*copy matrix: */
				matrix=xNew<bool>(M*N);
				memcpy(matrix,bmatrix,(M*N)*sizeof(bool));
			}

			else
				_error_("unrecognized bool pyarray type in input!");

			if (py_matrix2)
				delete(py_matrix2);
		}
		else
			matrix=NULL;
	}

	else {
		M=1;
		N=1;
		matrix=xNew<bool>(M*N);
		FetchData(&(matrix[0]),py_matrix);
	}

	/*output: */
	if(pM)*pM=M;
	if(pN)*pN=N;
	if(pmatrix)*pmatrix=matrix;
}
/*}}}*/
/*FUNCTION FetchData(double** pvector,int* pM, PyObject* py_vector){{{*/
void FetchData(double** pvector,int* pM,PyObject* py_vector){

	/*output: */
	double* dvector=NULL;
	double* vector=NULL;
	int M;
	int ndim;
	npy_intp*  dims=NULL;

	/*intermediary:*/
	long* lvector=NULL;
	bool* bvector=NULL;
	int i;
	PyObject* py_vector2=NULL;

	if     (PyArray_Check((PyArrayObject*)py_vector)) {
		/*retrieve dimensions: */
		ndim=PyArray_NDIM((const PyArrayObject*)py_vector);
		if      (ndim==1) {
			dims=PyArray_DIMS((PyArrayObject*)py_vector);
			M=dims[0]; 
		}
		else if (ndim==2) {
			dims=PyArray_DIMS((PyArrayObject*)py_vector);
			if (dims[1]==1)
				M=dims[0]; 
			else
				_error_("expecting an Mx1 matrix or M vector in input!");
		}
		else
			_error_("expecting an Mx1 matrix or M vector in input!");

		if (M) {
			if (!PyArray_ISCONTIGUOUS((PyArrayObject*)py_vector)) {
				py_vector2=PyArray_ContiguousFromAny(py_vector,NPY_DOUBLE,ndim,ndim);
				py_vector=py_vector2;
			}

			if      (PyArray_TYPE((PyArrayObject*)py_vector) == NPY_DOUBLE) {
				/*retrieve internal value: */
				dvector=(double*)PyArray_DATA((PyArrayObject*)py_vector);

				/*copy vector: */
				vector=xNew<double>(M);
				memcpy(vector,dvector,(M)*sizeof(double));
			}

			else if (PyArray_TYPE((PyArrayObject*)py_vector) == NPY_LONG) {
				/*retrieve internal value: */
				lvector=(long*)PyArray_DATA((PyArrayObject*)py_vector);

				/*transform into double vector: */
				vector=xNew<double>(M);
				for(i=0;i<M;i++)vector[i]=(double)lvector[i];
			}

			else if (PyArray_TYPE((PyArrayObject*)py_vector) == NPY_BOOL) {
				/*retrieve internal value: */
				bvector=(bool*)PyArray_DATA((PyArrayObject*)py_vector);

				/*transform into double vector: */
				vector=xNew<double>(M);
				for(i=0;i<M;i++)vector[i]=(double)bvector[i];
			}

			else
				_error_("unrecognized double pyarray type in input!");

			/* Causing a seg fault.
			if (py_vector2)
				delete(py_vector2);
			*/
		}
		else
			vector=NULL;
	}

	else {
		M=1;
		vector=xNew<double>(M);
		FetchData(&(vector[0]),py_vector);
	}

	/*output: */
	if(pM)*pM=M;
	if(pvector)*pvector=vector;
}
/*}}}*/
/*FUNCTION FetchData(int** pvector,int* pM, PyObject* py_vector){{{*/
void FetchData(int** pvector,int* pM,PyObject* py_vector){

	/*output: */
	int* vector=NULL;
	int M;
	int ndim;
	npy_intp*  dims=NULL;

	/*intermediary:*/
	long*   lvector=NULL;
	bool*   bvector=NULL;
	double* dvector=NULL;
	int i;
	PyObject* py_vector2=NULL;

	if     (PyArray_Check((PyArrayObject*)py_vector)) {
		/*retrieve dimensions: */
		ndim=PyArray_NDIM((const PyArrayObject*)py_vector);
		if      (ndim==1) {
			dims=PyArray_DIMS((PyArrayObject*)py_vector);
			M=dims[0]; 
		}
		else if (ndim==2) {
			dims=PyArray_DIMS((PyArrayObject*)py_vector);
			if (dims[1]==1)
				M=dims[0]; 
			else
				_error_("expecting an Mx1 matrix or M vector in input!");
		}
		else
			_error_("expecting an Mx1 matrix or M vector in input!");

		if (M) {
			if (!PyArray_ISCONTIGUOUS((PyArrayObject*)py_vector)) {
				py_vector2=PyArray_ContiguousFromAny(py_vector,NPY_LONG,ndim,ndim);
				py_vector=py_vector2;
			}

			if      (PyArray_TYPE((PyArrayObject*)py_vector) == NPY_DOUBLE) {
				/*retrieve internal value: */
				dvector=(double*)PyArray_DATA((PyArrayObject*)py_vector);

				/*transform into int vector: */
				vector=xNew<int>(M);
				for(i=0;i<M;i++)vector[i]=(int)lvector[i];
			}

			else if (PyArray_TYPE((PyArrayObject*)py_vector) == NPY_LONG) {
				/*retrieve internal value: */
				lvector=(long*)PyArray_DATA((PyArrayObject*)py_vector);

				/*transform into int vector: */
				vector=xNew<int>(M);
				for(i=0;i<M;i++)vector[i]=(int)lvector[i];
			}

			else if (PyArray_TYPE((PyArrayObject*)py_vector) == NPY_BOOL) {
				/*retrieve internal value: */
				bvector=(bool*)PyArray_DATA((PyArrayObject*)py_vector);

				/*transform into int vector: */
				vector=xNew<int>(M);
				for(i=0;i<M;i++)vector[i]=(int)bvector[i];
			}

			else
			 _error_("unrecognized int pyarray type in input!");

			if (py_vector2)
				delete(py_vector2);
		}
		else
			vector=NULL;
	}

	else {
		M=1;
		vector=xNew<int>(M);
		FetchData(&(vector[0]),py_vector);
	}

	/*output: */
	if(pM)*pM=M;
	if(pvector)*pvector=vector;
}
/*}}}*/
/*FUNCTION FetchData(bool** pvector,int* pM, PyObject* py_vector){{{*/
void FetchData(bool** pvector,int* pM,PyObject* py_vector){

	/*output: */
	bool* bvector=NULL;
	bool* vector=NULL;
	int M;
	int ndim;
	npy_intp*  dims=NULL;

	/*intermediary:*/
	double* dvector=NULL;
	long* lvector=NULL;
	int i;
	PyObject* py_vector2=NULL;

	if     (PyArray_Check((PyArrayObject*)py_vector)) {
		/*retrieve dimensions: */
		ndim=PyArray_NDIM((const PyArrayObject*)py_vector);
		if      (ndim==1) {
			dims=PyArray_DIMS((PyArrayObject*)py_vector);
			M=dims[0]; 
		}
		else if (ndim==2) {
			dims=PyArray_DIMS((PyArrayObject*)py_vector);
			if (dims[1]==1)
				M=dims[0]; 
			else
				_error_("expecting an Mx1 matrix or M vector in input!");
		}
		else
			_error_("expecting an Mx1 matrix or M vector in input!");

		if (M) {
			if (!PyArray_ISCONTIGUOUS((PyArrayObject*)py_vector)) {
				py_vector2=PyArray_ContiguousFromAny(py_vector,NPY_BOOL,ndim,ndim);
				py_vector=py_vector2;
			}

			if      (PyArray_TYPE((PyArrayObject*)py_vector) == NPY_DOUBLE) {
				/*retrieve internal value: */
				dvector=(double*)PyArray_DATA((PyArrayObject*)py_vector);

				/*transform into bool vector: */
				vector=xNew<bool>(M);
				for(i=0;i<M;i++)vector[i]=(bool)dvector[i];
			}

			else if (PyArray_TYPE((PyArrayObject*)py_vector) == NPY_LONG) {
				/*retrieve internal value: */
				lvector=(long*)PyArray_DATA((PyArrayObject*)py_vector);

				/*transform into bool vector: */
				vector=xNew<bool>(M);
				for(i=0;i<M;i++)vector[i]=(bool)lvector[i];
			}

			else if (PyArray_TYPE((PyArrayObject*)py_vector) == NPY_BOOL) {
				/*retrieve internal value: */
				bvector=(bool*)PyArray_DATA((PyArrayObject*)py_vector);

				/*copy vector: */
				vector=xNew<bool>(M);
				memcpy(vector,bvector,(M)*sizeof(bool));
			}

			else
				_error_("unrecognized bool pyarray type in input!");

			if (py_vector2)
				delete(py_vector2);
		}
		else
			vector=NULL;
	}

	else {
		M=1;
		vector=xNew<bool>(M);
		FetchData(&(vector[0]),py_vector);
	}

	/*output: */
	if(pM)*pM=M;
	if(pvector)*pvector=vector;
}
/*}}}*/

/*ISSM objects*/
/*FUNCTION FetchData(BamgGeom** pbamggeom,PyObject* py_dict){{{*/
void FetchData(BamgGeom** pbamggeom,PyObject* py_dict){

	/*Initialize output*/
	BamgGeom* bamggeom=new BamgGeom();

	/*Fetch all fields*/
	FetchData(&bamggeom->Vertices,&bamggeom->VerticesSize[0],&bamggeom->VerticesSize[1],PyDict_GetItemString(py_dict,"Vertices"));
	FetchData(&bamggeom->Edges, &bamggeom->EdgesSize[0], &bamggeom->EdgesSize[1], PyDict_GetItemString(py_dict,"Edges"));
	FetchData(&bamggeom->Corners, &bamggeom->CornersSize[0], &bamggeom->CornersSize[1], PyDict_GetItemString(py_dict,"Corners"));
	FetchData(&bamggeom->RequiredVertices,&bamggeom->RequiredVerticesSize[0],&bamggeom->RequiredVerticesSize[1],PyDict_GetItemString(py_dict,"RequiredVertices"));
	FetchData(&bamggeom->RequiredEdges, &bamggeom->RequiredEdgesSize[0], &bamggeom->RequiredEdgesSize[1], PyDict_GetItemString(py_dict,"RequiredEdges"));
	FetchData(&bamggeom->CrackedEdges,&bamggeom->CrackedEdgesSize[0],&bamggeom->CrackedEdgesSize[1],PyDict_GetItemString(py_dict,"CrackedEdges"));
	FetchData(&bamggeom->SubDomains,&bamggeom->SubDomainsSize[0],&bamggeom->SubDomainsSize[1],PyDict_GetItemString(py_dict,"SubDomains"));

	/*Assign output pointers:*/
	*pbamggeom=bamggeom;
}
/*}}}*/
/*FUNCTION FetchData(BamgMesh** pbamgmesh,PyObject* py_dict){{{*/
void FetchData(BamgMesh** pbamgmesh,PyObject* py_dict){

	/*Initialize output*/
	BamgMesh* bamgmesh=new BamgMesh();

	/*Fetch all fields*/
	FetchData(&bamgmesh->Vertices,&bamgmesh->VerticesSize[0],&bamgmesh->VerticesSize[1],PyDict_GetItemString(py_dict,"Vertices"));
	FetchData(&bamgmesh->Edges, &bamgmesh->EdgesSize[0], &bamgmesh->EdgesSize[1], PyDict_GetItemString(py_dict,"Edges"));
	FetchData(&bamgmesh->Triangles, &bamgmesh->TrianglesSize[0], &bamgmesh->TrianglesSize[1], PyDict_GetItemString(py_dict,"Triangles"));
	FetchData(&bamgmesh->CrackedEdges,&bamgmesh->CrackedEdgesSize[0],&bamgmesh->CrackedEdgesSize[1],PyDict_GetItemString(py_dict,"CrackedEdges"));
	FetchData(&bamgmesh->VerticesOnGeomEdge,&bamgmesh->VerticesOnGeomEdgeSize[0],&bamgmesh->VerticesOnGeomEdgeSize[1],PyDict_GetItemString(py_dict,"VerticesOnGeomEdge"));
	FetchData(&bamgmesh->VerticesOnGeomVertex,&bamgmesh->VerticesOnGeomVertexSize[0],&bamgmesh->VerticesOnGeomVertexSize[1],PyDict_GetItemString(py_dict,"VerticesOnGeomVertex"));
	FetchData(&bamgmesh->EdgesOnGeomEdge, &bamgmesh->EdgesOnGeomEdgeSize[0], &bamgmesh->EdgesOnGeomEdgeSize[1], PyDict_GetItemString(py_dict,"EdgesOnGeomEdge"));
	FetchData(&bamgmesh->IssmSegments,&bamgmesh->IssmSegmentsSize[0],&bamgmesh->IssmSegmentsSize[1],PyDict_GetItemString(py_dict,"IssmSegments"));

	/*Assign output pointers:*/
	*pbamgmesh=bamgmesh;
}
/*}}}*/
/*FUNCTION FetchData(BamgOpts** pbamgopts,PyObject* py_dict){{{*/
void FetchData(BamgOpts** pbamgopts,PyObject* py_dict){

	/*Initialize output*/
	BamgOpts* bamgopts=new BamgOpts();

	/*Fetch all fields*/
	FetchData(&bamgopts->anisomax,PyDict_GetItemString(py_dict,"anisomax"));
	FetchData(&bamgopts->cutoff,PyDict_GetItemString(py_dict,"cutoff"));
	FetchData(&bamgopts->coeff,PyDict_GetItemString(py_dict,"coeff"));
	FetchData(&bamgopts->errg,PyDict_GetItemString(py_dict,"errg"));
	FetchData(&bamgopts->gradation,PyDict_GetItemString(py_dict,"gradation"));
	FetchData(&bamgopts->Hessiantype,PyDict_GetItemString(py_dict,"Hessiantype"));
	FetchData(&bamgopts->MaxCornerAngle,PyDict_GetItemString(py_dict,"MaxCornerAngle"));
	FetchData(&bamgopts->maxnbv,PyDict_GetItemString(py_dict,"maxnbv"));
	FetchData(&bamgopts->maxsubdiv,PyDict_GetItemString(py_dict,"maxsubdiv"));
	FetchData(&bamgopts->Metrictype,PyDict_GetItemString(py_dict,"Metrictype"));
	FetchData(&bamgopts->nbjacobi,PyDict_GetItemString(py_dict,"nbjacobi"));
	FetchData(&bamgopts->nbsmooth,PyDict_GetItemString(py_dict,"nbsmooth"));
	FetchData(&bamgopts->omega,PyDict_GetItemString(py_dict,"omega"));
	FetchData(&bamgopts->power,PyDict_GetItemString(py_dict,"power"));
	FetchData(&bamgopts->random,PyDict_GetItemString(py_dict,"random"));
	FetchData(&bamgopts->verbose,PyDict_GetItemString(py_dict,"verbose"));

	FetchData(&bamgopts->Crack,PyDict_GetItemString(py_dict,"Crack"));
	FetchData(&bamgopts->geometricalmetric,PyDict_GetItemString(py_dict,"geometricalmetric"));
	FetchData(&bamgopts->KeepVertices,PyDict_GetItemString(py_dict,"KeepVertices"));
	FetchData(&bamgopts->splitcorners,PyDict_GetItemString(py_dict,"splitcorners"));

	FetchData(&bamgopts->hmin,PyDict_GetItemString(py_dict,"hmin"));
	FetchData(&bamgopts->hmax,PyDict_GetItemString(py_dict,"hmax"));
	FetchData(&bamgopts->hminVertices,&bamgopts->hminVerticesSize[0],&bamgopts->hminVerticesSize[1],PyDict_GetItemString(py_dict,"hminVertices"));
	FetchData(&bamgopts->hmaxVertices,&bamgopts->hmaxVerticesSize[0],&bamgopts->hmaxVerticesSize[1],PyDict_GetItemString(py_dict,"hmaxVertices"));
	FetchData(&bamgopts->hVertices,&bamgopts->hVerticesSize[0],&bamgopts->hVerticesSize[1],PyDict_GetItemString(py_dict,"hVertices"));
	FetchData(&bamgopts->metric,&bamgopts->metricSize[0],&bamgopts->metricSize[1],PyDict_GetItemString(py_dict,"metric"));
	FetchData(&bamgopts->field,&bamgopts->fieldSize[0],&bamgopts->fieldSize[1],PyDict_GetItemString(py_dict,"field"));
	FetchData(&bamgopts->err,&bamgopts->errSize[0],&bamgopts->errSize[1],PyDict_GetItemString(py_dict,"err"));

	/*Additional checks*/
	bamgopts->Check();

	/*Assign output pointers:*/
	*pbamgopts=bamgopts;
}
/*}}}*/
/*FUNCTION FetchData(Options** poptions,int istart, int nrhs,PyObject* py_tuple){{{*/
void FetchData(Options** poptions,int istart, int nrhs,PyObject* py_tuple){

	char   *name   = NULL;
	Option *option = NULL;

	/*Initialize output*/
	Options* options=new Options();

	/*Fetch all options*/
	for (int i=istart; i<nrhs; i=i+2){
		#if _PYTHON_MAJOR_ >= 3
		if (!PyUnicode_Check(PyTuple_GetItem(py_tuple,(Py_ssize_t)i))) _error_("Argument " << i+1 << " must be name of option");
		#else
		if (!PyString_Check(PyTuple_GetItem(py_tuple,(Py_ssize_t)i))) _error_("Argument " << i+1 << " must be name of option");
		#endif

		FetchData(&name,PyTuple_GetItem(py_tuple,(Py_ssize_t)i));
		if(i+1 == nrhs) _error_("Argument " << i+2 << " must exist and be value of option \"" << name << "\".");

		_printf0_("FetchData for Options not implemented yet, ignoring option \"" << name << "\"!\n");

//		option=(Option*)OptionParse(name,&PyTuple_GetItem(py_tuple,(Py_ssize_t)(i+1)));
//		options->AddOption(option);
//		option=NULL;
	}

	/*Assign output pointers:*/
	*poptions=options;
}
/*}}}*/
/*FUNCTION FetchData(Contours** pcontours,PyObject* py_list){{{*/
void FetchData(Contours** pcontours,PyObject* py_list){

	int              numcontours,test1,test2;
	char            *contourname = NULL;
	Contours         *contours    = NULL;
	Contour<double> *contouri    = NULL;
	PyObject        *py_dicti    = NULL;
	PyObject        *py_item     = NULL;

	#if _PYTHON_MAJOR_ >= 3
	if (PyUnicode_Check(py_list)){
		FetchData(&contourname,py_list);
		contours=ExpRead<double>(contourname);
	}
	#else
	if (PyString_Check(py_list)){
		FetchData(&contourname,py_list);
		contours=ExpRead<double>(contourname);
	}
	#endif
	else if(PyList_Check(py_list)){

		contours=new Contours();
		numcontours=(int)PyList_Size(py_list);

		for(int i=0;i<numcontours;i++){

			contouri=new Contour<double>();
			py_dicti=PyList_GetItem(py_list,(Py_ssize_t)i);

			py_item = PyDict_GetItemString(py_dicti,"nods");
			if(!py_item) _error_("input structure does not have a 'nods' field");
			FetchData(&contouri->nods,py_item);

			py_item = PyDict_GetItemString(py_dicti,"x");
			if(!py_item) _error_("input structure does not have a 'x' field");
			FetchData(&contouri->x,&test1,&test2,py_item);
			if(test1!=contouri->nods || test2!=1) _error_("field x should be of size ["<<contouri->nods<<" 1]");

			py_item = PyDict_GetItemString(py_dicti,"y");
			if(!py_item) _error_("input structure does not have a 'y' field");
			FetchData(&contouri->y,&test1,&test2,py_item);
			if(test1!=contouri->nods || test2!=1) _error_("field y should be of size ["<<contouri->nods<<" 1]");

			contours->AddObject(contouri);
		}
	}
	else{
		_error_("Contour is neither a string nor a structure and cannot be loaded");
	}

	/*clean-up and assign output pointer*/
	xDelete<char>(contourname);
	*pcontours=contours;
}
/*}}}*/

/*Python version dependent: */
#if _PYTHON_MAJOR_ >= 3 
/*FUNCTION FetchData(char** pstring,PyObject* py_unicode){{{*/
void FetchData(char** pstring,PyObject* py_unicode){

	PyObject* py_bytes;
	char* string=NULL;

	/*convert to bytes format: */
	PyUnicode_FSConverter(py_unicode,&py_bytes);

	/*convert from bytes to string: */
	string=PyBytes_AS_STRING(py_bytes);

	/*copy string (note strlen does not include trailing NULL): */
	*pstring=xNew<char>(strlen(string)+1);
	memcpy(*pstring,string,(strlen(string)+1)*sizeof(char));
}
/*}}}*/
#else
/*FUNCTION FetchData(char** pstring,PyObject* py_string){{{*/
void FetchData(char** pstring,PyObject* py_string){

	char* string=NULL;

	/*extract internal string: */
	string=PyString_AsString(py_string);

	/*copy string (note strlen does not include trailing NULL): */
	*pstring=xNew<char>(strlen(string)+1);
	memcpy(*pstring,string,(strlen(string)+1)*sizeof(char));
}
/*}}}*/
#endif
