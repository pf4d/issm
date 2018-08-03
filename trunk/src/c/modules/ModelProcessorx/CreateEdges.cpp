/*!\file:  CreateEdges.cpp
 * \brief: create edges from 2d mesh
 */ 

#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void CreateEdges(IoModel* iomodel){/*{{{*/

	/*If edges are already present, exit*/
	if(iomodel->edges) return;

	/*Check Iomodel properties*/
	if(iomodel->numberofvertices<3) _error_("not enough elements in mesh");
	_assert_(iomodel->elements);

	/*Intermediaries*/
	bool exist;
	int  i,j,v1,v2,v3;
	int  maxnbe,nbe,elementnbe,elementnbv;
	int *elementedges         = NULL;
	int *elementedges_markers = NULL;

	/*Mesh dependent variables*/
	switch(iomodel->meshelementtype){
		case TriaEnum:
			elementnbv = 3;
			elementnbe = 3;
			elementedges         = xNew<int>(elementnbe*2);
			elementedges_markers = xNew<int>(elementnbe);
			elementedges[2*0+0] = 1; elementedges[2*0+1] = 2; elementedges_markers[0] = 1;
			elementedges[2*1+0] = 2; elementedges[2*1+1] = 0; elementedges_markers[1] = 1;
			elementedges[2*2+0] = 0; elementedges[2*2+1] = 1; elementedges_markers[2] = 1;
			break;
		case TetraEnum:
			elementnbv = 4;
			elementnbe = 6;
			elementedges         = xNew<int>(elementnbe*2);
			elementedges_markers = xNew<int>(elementnbe);
			elementedges[2*0+0] = 1; elementedges[2*0+1] = 2; elementedges_markers[0] = 1;
			elementedges[2*1+0] = 0; elementedges[2*1+1] = 2; elementedges_markers[1] = 1;
			elementedges[2*2+0] = 0; elementedges[2*2+1] = 1; elementedges_markers[2] = 1;
			elementedges[2*3+0] = 1; elementedges[2*3+1] = 3; elementedges_markers[3] = 1;
			elementedges[2*4+0] = 2; elementedges[2*4+1] = 3; elementedges_markers[4] = 1;
			elementedges[2*5+0] = 0; elementedges[2*5+1] = 3; elementedges_markers[5] = 1;
			break;
		case PentaEnum:
			elementnbv = 6;
			elementnbe = 9;
			elementedges         = xNew<int>(elementnbe*2);
			elementedges_markers = xNew<int>(elementnbe);
			elementedges[2*0+0] = 0; elementedges[2*0+1] = 3; elementedges_markers[0] = 2;
			elementedges[2*1+0] = 1; elementedges[2*1+1] = 4; elementedges_markers[1] = 2;
			elementedges[2*2+0] = 2; elementedges[2*2+1] = 5; elementedges_markers[2] = 2;
			elementedges[2*3+0] = 1; elementedges[2*3+1] = 2; elementedges_markers[3] = 1;
			elementedges[2*4+0] = 2; elementedges[2*4+1] = 0; elementedges_markers[4] = 1;
			elementedges[2*5+0] = 0; elementedges[2*5+1] = 1; elementedges_markers[5] = 1;
			elementedges[2*6+0] = 4; elementedges[2*6+1] = 5; elementedges_markers[6] = 1;
			elementedges[2*7+0] = 5; elementedges[2*7+1] = 3; elementedges_markers[7] = 1;
			elementedges[2*8+0] = 3; elementedges[2*8+1] = 4; elementedges_markers[8] = 1;
			break;
		default:
		_error_("mesh dimension not supported yet");
	}

	/*Maximum number of edges*/
	maxnbe = elementnbe*iomodel->numberofelements;

	/*Initialize intermediaries*/
	int *edgestemp                 = xNew<int>(maxnbe*3);                             /*format: [vertex1 vertex2 marker]       */
	int *element_edge_connectivity = xNew<int>(iomodel->numberofelements*elementnbe); /*format: [edge1 edge2 ... edgen] */

	/*Initialize chain*/
	int* head_minv = xNew<int>(iomodel->numberofvertices);
	int* next_edge = xNew<int>(maxnbe);
	for(i=0;i<iomodel->numberofvertices;i++) head_minv[i]=-1;

	/*Initialize number of edges*/
	nbe = 0;

	for(i=0;i<iomodel->numberofelements;i++){
		for(j=0;j<elementnbe;j++){

			/*Get the two indices of the edge number j of the ith element*/
			v1 = iomodel->elements[i*elementnbv+elementedges[2*j+0]]-1; _assert_(v1>=0 & v1<iomodel->numberofvertices);
			v2 = iomodel->elements[i*elementnbv+elementedges[2*j+1]]-1; _assert_(v2>=0 & v2<iomodel->numberofvertices);

			/*v1 and v2 must be sorted*/
			if(v2<v1){
				v3=v2; v2=v1; v1=v3;
			}

			/*This edge a priori has not been processed yet*/
			exist = false;

			/*Go through all processed edges connected to v1 and check whether we have seen this edge yet*/
			for(int e=head_minv[v1]; e!=-1; e=next_edge[e]){
				if(edgestemp[e*3+1]==v2+1){
					exist = true;
					element_edge_connectivity[i*elementnbe+j]=e;
					break;
				}
			}

			/*If this edge is new, add it to the lists*/
			if(!exist){
				_assert_(nbe<maxnbe);

				/*Update edges*/
				edgestemp[nbe*3+0] = v1+1;
				edgestemp[nbe*3+1] = v2+1;
				edgestemp[nbe*3+2] = elementedges_markers[j];

				/*Update Connectivity*/
				element_edge_connectivity[i*elementnbe+j]=nbe;

				/*Update chain*/
				next_edge[nbe] = head_minv[v1];
				head_minv[v1]  = nbe;

				/*Increase number of edges*/
				nbe++;
			}
		}
	}

	/*Clean up*/
	xDelete<int>(head_minv);
	xDelete<int>(next_edge);
	xDelete<int>(elementedges_markers);

	/*Create final edges*/
	int* edges = xNew<int>(nbe*3); /*format: [vertex1 vertex2 marker]*/
	for(int i=0;i<3*nbe;i++) edges[i] = edgestemp[i];

	/*Clean up*/
	xDelete<int>(edgestemp);
	xDelete<int>(elementedges);

	/*Assign output pointers*/
	iomodel->edges                     = edges;
	iomodel->elementtoedgeconnectivity = element_edge_connectivity;
	iomodel->numberofedges             = nbe;
}/*}}}*/
void EdgeOnBoundaryFlags(bool** pflags,IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	bool isv1,isv2;
	int  facenbv,v1,v2;
	int  id_edge,id_element;
	int  elementnbe;

	/*Mesh dependent variables*/
	switch(iomodel->meshelementtype){
		case TriaEnum:  elementnbe = 3; break;
		case TetraEnum: elementnbe = 6; break;
		case PentaEnum: elementnbe = 9; break;
		default:        _error_("mesh dimension not supported yet");
	}

	/*Get edges and allocate output*/
	if(!iomodel->edges) CreateEdges(iomodel);
	bool* flags = xNewZeroInit<bool>(iomodel->numberofedges);

	if(iomodel->domaindim==2){

		/*Count how many times an edge is found in elementtoedgeconnectivity*/
		int* counter = xNewZeroInit<int>(iomodel->numberofedges);
		for(int i=0;i<iomodel->numberofelements;i++){
			for(int j=0;j<elementnbe;j++){
				counter[iomodel->elementtoedgeconnectivity[elementnbe*i+j]] += 1;
			}
		}

		/*Now, loop over the egdes, whenever it is not connected to a second element, the edge is on boundary*/
		for(int i=0;i<iomodel->numberofedges;i++){
			if(counter[i]==1) flags[i]=true;
		}

		/*Clean up*/
		xDelete<int>(counter);
	}
	else if(iomodel->domaindim==3){

		/*Get faces*/
		if(!iomodel->faces) CreateFaces(iomodel);

		/*Now, loop over the faces, whenever it is not connected to a second element, all edges are on boundary*/
		for(int id_face=0;id_face<iomodel->numberoffaces;id_face++){

			if(iomodel->faces[id_face*iomodel->facescols+1]==-1){

				/*The face is connected to the element e only*/
				id_element = iomodel->faces[id_face*iomodel->facescols+0]-1;
				facenbv    = iomodel->faces[id_face*iomodel->facescols+3];

				/*Get all edges for this element*/
				for(int edge = 0; edge<elementnbe; edge++){

					id_edge     = iomodel->elementtoedgeconnectivity[elementnbe*id_element+edge];
					v1          = iomodel->edges[id_edge*3+0];
					v2          = iomodel->edges[id_edge*3+1];

					/*Test if v1 is in the face*/
					isv1=false;
					for(int i=0;i<facenbv;i++){
						if(iomodel->faces[id_face*iomodel->facescols+4+i] == v1){
							isv1 = true; break;
						}
					}
					if(!isv1) continue;

					/*test if v2 is in the face*/
					isv2=false;
					for(int i=0;i<facenbv;i++){
						if(iomodel->faces[id_face*iomodel->facescols+4+i] == v2){
							isv2 = true; break;
						}
					}

					/*If v1 and v2 are found, this edge is on boundary*/
					if(isv2) flags[id_edge] = true;
				}
			}
		}
	}
	else{
		_error_("dimension not supported");
	}

	/*Clean up and return*/
	*pflags = flags;
}/*}}}*/
