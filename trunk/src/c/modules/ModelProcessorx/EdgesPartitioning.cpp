/*!\file:  EdgesPartitioning.cpp
 * \brief: partition elements and nodes and vertices
 */ 

#include <string.h>
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void EdgesPartitioning(bool** pmy_edges,IoModel* iomodel){

	/*Intermediaries*/
	int elementnbe;

	/*Get edges and elements*/
	CreateEdges(iomodel);
	_assert_(iomodel->elementtoedgeconnectivity);

	/*Mesh dependent variables*/
	switch(iomodel->meshelementtype){
		case TriaEnum:  elementnbe = 3; break;
		case TetraEnum: elementnbe = 6; break;
		case PentaEnum: elementnbe = 9; break;
		default: _error_("mesh dimension not supported yet");
	}

	/*output: */
	bool* my_edges=xNewZeroInit<bool>(iomodel->numberofedges);

	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			for(int j=0;j<elementnbe;j++){
				my_edges[iomodel->elementtoedgeconnectivity[i*elementnbe+j]] = true;
			}
		}
	}

	/*Free data and assign output pointers */
	*pmy_edges=my_edges;
}
