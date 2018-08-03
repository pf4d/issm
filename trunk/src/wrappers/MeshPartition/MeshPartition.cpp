/*!\file:  MeshPartition.cpp
 * \brief: partition mesh according to number of areas, using Metis library.
*/

#include "./MeshPartition.h"

void MeshPartitionUsage(void){/*{{{*/
	_printf_("   usage:\n");
	_printf_("   [element_partitioning,node_partitioning]=MeshPartition(md.mesh,numpartitions)");
	_printf_("   where:\n");
	_printf_("      element_partitioning is a vector of partitioning area numbers, for every element.\n");
	_printf_("      node_partitioning is a vector of partitioning area numbers, for every node.\n");
	_printf_("\n");
}/*}}}*/
WRAPPER(MeshPartition_python){

	/*Indexing: */
	int i,j;

	/* required input: */
	int  meshelementtype;
	int  numberofelements;
	int  numberofvertices;
	int *elements         = NULL;
	int  elements_width;

	int numberofelements2d;
	int numberofvertices2d;
	int* elements2d=NULL;

	int numberoflayers;
	int numareas=1;

	/* output: */
	int    *int_element_partitioning = NULL;
	int    *int_node_partitioning    = NULL;
	double *element_partitioning     = NULL;
	double *node_partitioning        = NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	CheckNumMatlabArguments(nlhs,NLHS,nrhs,NRHS,__FUNCT__,&MeshPartitionUsage);

	/*Fetch data: */
	FetchData(&numberofelements,mxGetAssignedField(MESH,0,"numberofelements"));
	FetchData(&numberofvertices,mxGetAssignedField(MESH,0,"numberofvertices"));
	FetchData(&elements,NULL,&elements_width,mxGetAssignedField(MESH,0,"elements"));

	if(strcmp(mxGetClassName(MESH),"mesh3dprisms")==0){
		meshelementtype = PentaEnum;
		FetchData(&numberofelements2d,mxGetAssignedField(MESH,0,"numberofelements2d"));
		FetchData(&numberofvertices2d,mxGetAssignedField(MESH,0,"numberofvertices2d"));
		FetchData(&elements2d,NULL,NULL,mxGetAssignedField(MESH,0,"elements2d"));
		FetchData(&numberoflayers,mxGetAssignedField(MESH,0,"numberoflayers"));
	}
	else if(strcmp(mxGetClassName(MESH),"mesh2dhorizontal")==0){
		meshelementtype = TriaEnum;
		numberoflayers=1;
	}
	else if(strcmp(mxGetClassName(MESH),"mesh2dvertical")==0){
		meshelementtype = TriaEnum;
		numberoflayers=1;
	}
	else{
		_error_("Mesh type "<<mxGetClassName(MESH)<<" not supported yet");
	}
	FetchData(&numareas,NUMAREAS);

	/*Run partitioning algorithm based on a "clever" use of the Metis partitioner: */
	MeshPartitionx(&int_element_partitioning,&int_node_partitioning,numberofelements,numberofvertices,elements,
		numberofelements2d,numberofvertices2d,elements2d,numberoflayers,elements_width,meshelementtype,numareas);

	/*Post process node_partitioning and element_partitioning to be in double format. Metis needed them in int* format: */
	element_partitioning=xNew<double>(numberofelements);
	for (i=0;i<numberofelements;i++){
		element_partitioning[i]=(double)int_element_partitioning[i]+1; //Metis indexing from 0, matlab from 1.
	}

	node_partitioning=xNew<double>(numberofvertices);
	for (i=0;i<numberofvertices;i++){
		node_partitioning[i]=(double)int_node_partitioning[i]+1; //Metis indexing from 0, matlab from 1.
	}

	/*Write data:*/
	WriteData(ELEMENTPARTITIONING,element_partitioning,numberofelements);
	WriteData(NODEPARTITIONING,node_partitioning,numberofvertices);

	/*Free ressources:*/
	xDelete<int>(elements);
	xDelete<int>( elements2d);
	xDelete<int>(int_element_partitioning);
	xDelete<int>(int_node_partitioning);
	xDelete<double>(element_partitioning);
	xDelete<double>(node_partitioning);

	/*end module: */
	MODULEEND();
}
