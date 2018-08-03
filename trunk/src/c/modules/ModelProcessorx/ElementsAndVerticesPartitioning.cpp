/*!\file:  ElementsAndVerticesPartitioning.cpp
 * \brief: partition elements and nodes and vertices
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <string.h>
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "../MeshPartitionx/MeshPartitionx.h"
#include "../ModelProcessorx/ModelProcessorx.h"

void  ElementsAndVerticesPartitioning(bool** pmy_elements, int** pmy_vertices, IoModel* iomodel){

	int i,j;

	const int RIFTINFOSIZE = 12;
	int numberofelements2d;
	int numberofvertices2d;
	int numlayers;
	int numrifts;
	int numvertex_pairing;

	/*output: */
	bool *my_elements = NULL;
	int  *my_vertices = NULL;

	/*intermediary: */
	int        *epart          = NULL; //element partitioning.
	int        *npart          = NULL; //node partitioning.
	int         elements_width;        //number of columns in elements (2d->3, 3d->6)
	int         el1,el2;
	int        *elements2d     = NULL;
	int        *vertex_pairing = NULL;
	IssmDouble *riftinfo       = NULL;

	/*Get my_rank:*/
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*Fetch parameters: */
	iomodel->FindConstant(&numrifts,"md.rifts.numrifts");

	/*First, check that partitioning has not yet been carryed out. Just check whether my_elements pointers is not already assigned a value: */
	if(*pmy_elements)return;

	/*Number of vertices per elements, needed to correctly retrieve data: */
	/*Determine parallel partitioning of elements: we use Metis for now. First load the data, then partition*/
	switch(iomodel->meshelementtype){
		case TriaEnum:
			elements_width=3;
			numberofelements2d = 0;
			numberofvertices2d = 0;
			numlayers          = 0;
			break;
		case TetraEnum:
			elements_width=4;
			numberofelements2d = 0;
			numberofvertices2d = 0;
			numlayers          = 0;
			break;
		case PentaEnum:
			elements_width=6;
			iomodel->FetchData(&elements2d,NULL,NULL,"md.mesh.elements2d");
			iomodel->FindConstant(&numberofelements2d,"md.mesh.numberofelements2d");
			iomodel->FindConstant(&numberofvertices2d,"md.mesh.numberofvertices2d");
			iomodel->FindConstant(&numlayers,"md.mesh.numberoflayers");
			break;
		default:
			_error_("mesh elements "<< EnumToStringx(iomodel->meshelementtype) <<" not supported yet");
	}

	MeshPartitionx(&epart,&npart,iomodel->numberofelements,iomodel->numberofvertices,iomodel->elements,numberofelements2d,numberofvertices2d,elements2d,numlayers,elements_width,iomodel->meshelementtype,num_procs);

	/*Free elements2d: */
	xDelete<int>(elements2d);

	/*Deal with rifts, they have to be included into one partition only, not several: */
	if(numrifts){
		iomodel->FetchData(&riftinfo,&numrifts,NULL,"md.rifts.riftstruct");
		for(i=0;i<numrifts;i++){
			el1=reCast<int>(*(riftinfo+RIFTINFOSIZE*i+2))-1; //matlab indexing to c indexing
			el2=reCast<int>(*(riftinfo+RIFTINFOSIZE*i+3))-1; //matlab indexing to c indexing
			epart[el2]=epart[el1]; //ensures that this pair of elements will be in the same partition, as well as the corresponding vertices;
		}
		iomodel->DeleteData(riftinfo,"md.rifts.riftstruct");
	}

	/*Used later on: */
	my_vertices=xNewZeroInit<int>(iomodel->numberofvertices);
	my_elements=xNewZeroInit<bool>(iomodel->numberofelements);

	/*Start figuring out, out of the partition, which elements belong to this cpu: */
	for(i=0;i<iomodel->numberofelements;i++){

		/*!All elements have been partitioned above, only deal with elements for this cpu: */
		if(my_rank==epart[i]){ 
			my_elements[i]=true;
			/*Now that we are here, we can also start building the list of vertices belonging to this cpu partition: we use 
			 *the  element index to do this. For each element n, we know index[n][0:2] holds the indices (matlab indexing) 
			 into the vertices coordinates. If we start plugging 1 into my_vertices for each index[n][i] (i=0:2), then my_vertices 
			 will hold which vertices belong to this partition*/
			for(j=0;j<elements_width;j++){
				my_vertices[iomodel->elements[elements_width*i+j]-1]=1;
			}
		}
	}

	/*We might have vertex_pairing in which case, some vertices have to be cloned:
	 * penpair has 2 nodes that are poointing toward 2 vertices.
	 * The 2 vertices must be in the same cpu as the penpair*/
	iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.stressbalance.vertex_pairing");
	for(i=0;i<numvertex_pairing;i++){
		if(my_vertices[vertex_pairing[2*i+0]-1] && !my_vertices[vertex_pairing[2*i+1]-1]){
			my_vertices[vertex_pairing[2*i+1]-1]=2; //to know that these elements are not on the partition
		}
	}
	xDelete<int>(vertex_pairing);
	iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.masstransport.vertex_pairing");
	for(i=0;i<numvertex_pairing;i++){
		if(my_vertices[vertex_pairing[2*i+0]-1] && !my_vertices[vertex_pairing[2*i+1]-1]){
			my_vertices[vertex_pairing[2*i+1]-1]=2; //to know that these elements are not on the partition
		}
	}
	xDelete<int>(vertex_pairing);

	/*Free ressources:*/
	xDelete<int>(npart);
	xDelete<int>(epart);

	/*Assign output pointers:*/
	*pmy_elements=my_elements;
	*pmy_vertices=my_vertices;
}
