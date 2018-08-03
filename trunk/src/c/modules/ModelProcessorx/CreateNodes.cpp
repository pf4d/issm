/*
 * CreateNodes.c:
 */

#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "../MeshPartitionx/MeshPartitionx.h"
#include "./ModelProcessorx.h"

void CreateNodes(Nodes* nodes, IoModel* iomodel,int analysis,int finite_element,int approximation){

	/*Intermediaries*/
	int   i,j,counter,vnodes,lid=0;
	int   numberoffaces,elementnbv;
	int   id0 = iomodel->nodecounter;
	bool *my_faces = NULL;
	bool *my_edges = NULL;
	bool *my_nodes = NULL;
	Node *node     = NULL;

	switch(finite_element){
		case P1Enum:
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,approximation));
				}
			}
			break;

		case P1DGEnum:
			NodesPartitioning(&my_nodes,iomodel->my_elements,iomodel->my_vertices,iomodel,false);
			for(i=0;i<iomodel->numberofelements;i++){
				for(j=0;j<3;j++){
					if(my_nodes[3*i+j]){ 
						nodes->AddObject(new Node(id0+3*i+j+1,id0+3*i+j,lid++,iomodel->elements[+3*i+j]-1,iomodel,analysis,approximation));

					}
				}
			}
			break;

		case P1bubbleEnum:
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,approximation));
				}
			}
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					nodes->AddObject(new Node(id0+iomodel->numberofvertices+i+1,iomodel->numberofvertices+i,lid++,0,iomodel,analysis,approximation));
				}
			}
			break;

		case P1bubblecondensedEnum:
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,approximation));
				}
			}
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					node = new Node(id0+iomodel->numberofvertices+i+1,iomodel->numberofvertices+i,lid++,0,iomodel,analysis,approximation);
					node->HardDeactivate();
					nodes->AddObject(node);
				}
			}
			break;

		case P1xP2Enum:
			EdgesPartitioning(&my_edges,iomodel);
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,approximation));
				}
			}

			counter = iomodel->numberofvertices;
			for(i=0;i<iomodel->numberofedges;i++){
				if(iomodel->edges[i*3+2]==2){
					if(my_edges[i]){
						node = new Node(id0+iomodel->numberofvertices+i+1,counter+1,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
					}
					counter++;
				}
			}
			break;

		case P1xP3Enum:
			EdgesPartitioning(&my_edges,iomodel);
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,approximation));
				}
			}

			counter = iomodel->numberofvertices;
			for(i=0;i<iomodel->numberofedges;i++){
				if(iomodel->edges[i*3+2]==2){
					if(my_edges[i]){
						node = new Node(id0+iomodel->numberofvertices+2*i+1,counter+1,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
						node = new Node(id0+iomodel->numberofvertices+2*i+2,counter+2,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
					}
					counter=counter+2;
				}
			}
			break;

		case P2xP1Enum:
			EdgesPartitioning(&my_edges,iomodel);
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,approximation));
				}
			}

			counter = iomodel->numberofvertices;
			for(i=0;i<iomodel->numberofedges;i++){
				if(iomodel->edges[i*3+2]!=2){
					if(my_edges[i]){
						node = new Node(id0+iomodel->numberofvertices+i+1,counter+1,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
					}
					counter++;
				}
			}
			break;

		case P2Enum:
			EdgesPartitioning(&my_edges,iomodel);
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,approximation));
				}
			}
			for(i=0;i<iomodel->numberofedges;i++){
				if(my_edges[i]){
					nodes->AddObject(new Node(id0+iomodel->numberofvertices+i+1,iomodel->numberofvertices+i,lid++,0,iomodel,analysis,approximation));
				}
			}
			id0 = id0+iomodel->numberofvertices+iomodel->numberofedges;
	      if(iomodel->meshelementtype==PentaEnum){
				FacesPartitioning(&my_faces,iomodel);
				for(i=0;i<iomodel->numberoffaces;i++){
					if(iomodel->faces[i*iomodel->facescols+2]==2){/*Vertical quads*/
						if(my_faces[i]){
							node = new Node(id0+i+1,iomodel->numberofvertices+iomodel->numberofedges+i,lid++,0,iomodel,analysis,approximation);
							nodes->AddObject(node);
						}
					}
					else if(iomodel->faces[i*iomodel->facescols+2]==1){/*Triangular base/top*/
						/*Nothing*/
					}
					else{
						_error_("not supported");
					}
				}
			}
			break;
		case P2bubbleEnum:
			EdgesPartitioning(&my_edges,iomodel);
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,approximation));
				}
			}
			for(i=0;i<iomodel->numberofedges;i++){
				if(my_edges[i]){
					nodes->AddObject(new Node(id0+iomodel->numberofvertices+i+1,iomodel->numberofvertices+i,lid++,0,iomodel,analysis,approximation));
				}
			}
			id0 = id0+iomodel->numberofvertices+iomodel->numberofedges;
			if(iomodel->meshelementtype==PentaEnum){
				FacesPartitioning(&my_faces,iomodel);
				for(i=0;i<iomodel->numberoffaces;i++){
					if(iomodel->faces[i*iomodel->facescols+2]==2){/*Vertical quads*/
						if(my_faces[i]){
							node = new Node(id0+i+1,iomodel->numberofvertices+iomodel->numberofedges+i,lid++,0,iomodel,analysis,approximation);
							nodes->AddObject(node);
						}
					}
					else if(iomodel->faces[i*iomodel->facescols+2]==1){/*Triangular base/top*/
						/*Nothing*/
					}
					else{
						_error_("not supported");
					}
				}
				id0 = id0+iomodel->numberoffaces;
			}
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					nodes->AddObject(new Node(id0+i+1,id0-iomodel->nodecounter+i,lid++,0,iomodel,analysis,approximation));
				}
			}
			break;
		case P2xP4Enum:
			EdgesPartitioning(&my_edges,iomodel);
			FacesPartitioning(&my_faces,iomodel);
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,approximation));
				}
			}
			counter = iomodel->numberofvertices;
			for(i=0;i<iomodel->numberofedges;i++){
				if(iomodel->edges[i*3+2]==2){/*Vertical edges*/
					if(my_edges[i]){
						node = new Node(id0+iomodel->numberofvertices+3*i+1,counter+1,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
						node = new Node(id0+iomodel->numberofvertices+3*i+2,counter+2,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
						node = new Node(id0+iomodel->numberofvertices+3*i+3,counter+3,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
					}
					counter=counter+3;
				}
				else if(iomodel->edges[i*3+2]==1){/*Horizontal edges*/
					if(my_edges[i]){
						node = new Node(id0+iomodel->numberofvertices+3*i+1,counter+1,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
					}
					counter=counter+1;
				}
				else{
					_error_("not supported");
				}
			}
			id0 = id0+iomodel->numberofvertices+3*iomodel->numberofedges;
			for(i=0;i<iomodel->numberoffaces;i++){
				if(iomodel->faces[i*iomodel->facescols+2]==2){/*Vertical quads*/
					if(my_faces[i]){
						node = new Node(id0+3*i+1,counter+1,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
						node = new Node(id0+3*i+2,counter+2,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
						node = new Node(id0+3*i+3,counter+3,lid++,0,iomodel,analysis,approximation);
						nodes->AddObject(node);
					}
					counter=counter+3;
				}
				else if(iomodel->faces[i*iomodel->facescols+2]==1){/*Triangular base/top*/
					/*Nothing*/
				}
				else{
					_error_("not supported");
				}
			}
			break;

		/*Stokes elements*/
		case P1P1Enum:
			_assert_(approximation==FSApproximationEnum);
			/*P1 velocity*/
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,FSvelocityEnum));
				}
			}
			/*P1 pressure*/
			vnodes = id0+iomodel->numberofvertices;
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(vnodes+i+1,iomodel->numberofvertices+i,lid++,i,iomodel,analysis,FSpressureEnum));
				}
			}
			break;
		case P1P1GLSEnum:
			_assert_(approximation==FSApproximationEnum);
			/*P1 velocity*/
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,FSvelocityEnum));
				}
			}
			/*P1 pressure*/
			vnodes = id0+iomodel->numberofvertices;
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(vnodes+i+1,iomodel->numberofvertices+i,lid++,i,iomodel,analysis,FSpressureEnum));
				}
			}
			break;
		case MINIcondensedEnum:
			_assert_(approximation==FSApproximationEnum);
			/*P1+ velocity (bubble statically condensed)*/
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,FSvelocityEnum));
				}
			}
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					node = new Node(id0+iomodel->numberofvertices+i+1,iomodel->numberofvertices+i,lid++,0,iomodel,analysis,FSvelocityEnum);
					node->HardDeactivate();
					nodes->AddObject(node);
				}
			}
			/*P1 pressure*/
			vnodes = id0+iomodel->numberofvertices+iomodel->numberofelements;
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(vnodes+i+1,iomodel->numberofvertices+iomodel->numberofelements+i,lid++,i,iomodel,analysis,FSpressureEnum));
				}
			}
			break;
		case MINIEnum:
			_assert_(approximation==FSApproximationEnum);
			/*P1+ velocity*/
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,FSvelocityEnum));
				}
			}
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					nodes->AddObject(new Node(id0+iomodel->numberofvertices+i+1,iomodel->numberofvertices+i,lid++,0,iomodel,analysis,FSvelocityEnum));
				}
			}
			/*P1 pressure*/
			vnodes = id0+iomodel->numberofvertices+iomodel->numberofelements;
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(vnodes+i+1,iomodel->numberofvertices+iomodel->numberofelements+i,lid++,i,iomodel,analysis,FSpressureEnum));
				}
			}
			break;
		case TaylorHoodEnum:
		case XTaylorHoodEnum:
			_assert_(approximation==FSApproximationEnum);
			/*P2 velocity*/
			EdgesPartitioning(&my_edges,iomodel);
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,FSvelocityEnum));
				}
			}
			for(i=0;i<iomodel->numberofedges;i++){
				if(my_edges[i]){
					nodes->AddObject(new Node(id0+iomodel->numberofvertices+i+1,iomodel->numberofvertices+i,lid++,0,iomodel,analysis,FSvelocityEnum));
				}
			}
			id0 = id0+iomodel->numberofvertices+iomodel->numberofedges;
	      if(iomodel->meshelementtype==PentaEnum){
				FacesPartitioning(&my_faces,iomodel);
				for(i=0;i<iomodel->numberoffaces;i++){
					if(iomodel->faces[i*iomodel->facescols+2]==2){/*Vertical quads*/
						if(my_faces[i]){
							node = new Node(id0+i+1,iomodel->numberofvertices+iomodel->numberofedges+i,lid++,0,iomodel,analysis,FSvelocityEnum);
							nodes->AddObject(node);
						}
					}
					else if(iomodel->faces[i*iomodel->facescols+2]==1){/*Triangular base/top*/
						/*Nothing*/
					}
					else{
						_error_("not supported");
					}
				}
			}

			/*P1 pressure*/
	      if(iomodel->meshelementtype==PentaEnum){
				numberoffaces=iomodel->numberoffaces;
			}
			else{
				numberoffaces=0;
			}
			vnodes = id0+numberoffaces;
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(vnodes+i+1,iomodel->numberofvertices+iomodel->numberofedges+numberoffaces+i,lid++,i,iomodel,analysis,FSpressureEnum));
				}
			}
			break;
		case LATaylorHoodEnum:
			_assert_(approximation==FSApproximationEnum);
			/*P2 velocity*/
			EdgesPartitioning(&my_edges,iomodel);
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,FSvelocityEnum));
				}
			}
			for(i=0;i<iomodel->numberofedges;i++){
				if(my_edges[i]){
					nodes->AddObject(new Node(id0+iomodel->numberofvertices+i+1,iomodel->numberofvertices+i,lid++,0,iomodel,analysis,FSvelocityEnum));
				}
			}
			id0 = id0+iomodel->numberofvertices+iomodel->numberofedges;
	      if(iomodel->meshelementtype==PentaEnum){
				FacesPartitioning(&my_faces,iomodel);
				for(i=0;i<iomodel->numberoffaces;i++){
					if(iomodel->faces[i*iomodel->facescols+2]==2){/*Vertical quads*/
						if(my_faces[i]){
							node = new Node(id0+i+1,iomodel->numberofvertices+iomodel->numberofedges+i,lid++,0,iomodel,analysis,FSvelocityEnum);
							nodes->AddObject(node);
						}
					}
					else if(iomodel->faces[i*iomodel->facescols+2]==1){/*Triangular base/top*/
						/*Nothing*/
					}
					else{
						_error_("not supported");
					}
				}
			}

			/*No pressure*/
			break;
		case OneLayerP4zEnum:
			_assert_(approximation==FSApproximationEnum);
			EdgesPartitioning(&my_edges,iomodel);
			FacesPartitioning(&my_faces,iomodel);
			/*P2xP4 velocity*/
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,FSvelocityEnum));
				}
			}
			counter = iomodel->numberofvertices;
			for(i=0;i<iomodel->numberofedges;i++){
				if(iomodel->edges[i*3+2]==2){/*Vertical edges*/
					if(my_edges[i]){
						node = new Node(id0+iomodel->numberofvertices+3*i+1,counter+1,lid++,0,iomodel,analysis,FSvelocityEnum);
						nodes->AddObject(node);
						node = new Node(id0+iomodel->numberofvertices+3*i+2,counter+2,lid++,0,iomodel,analysis,FSvelocityEnum);
						nodes->AddObject(node);
						node = new Node(id0+iomodel->numberofvertices+3*i+3,counter+3,lid++,0,iomodel,analysis,FSvelocityEnum);
						nodes->AddObject(node);
					}
					counter=counter+3;
				}
				else if(iomodel->edges[i*3+2]==1){/*Horizontal edges*/
					if(my_edges[i]){
						node = new Node(id0+iomodel->numberofvertices+3*i+1,counter+1,lid++,0,iomodel,analysis,FSvelocityEnum);
						nodes->AddObject(node);
					}
					counter=counter+1;
				}
				else{
					_error_("not supported");
				}
			}
			id0 = id0+iomodel->numberofvertices+3*iomodel->numberofedges;
			for(i=0;i<iomodel->numberoffaces;i++){
				if(iomodel->faces[i*iomodel->facescols+2]==2){/*Vertical quads*/
					if(my_faces[i]){
						node = new Node(id0+3*i+1,counter+1,lid++,0,iomodel,analysis,FSvelocityEnum);
						nodes->AddObject(node);
						node = new Node(id0+3*i+2,counter+2,lid++,0,iomodel,analysis,FSvelocityEnum);
						nodes->AddObject(node);
						node = new Node(id0+3*i+3,counter+3,lid++,0,iomodel,analysis,FSvelocityEnum);
						nodes->AddObject(node);
					}
					counter=counter+3;
				}
				else if(iomodel->faces[i*iomodel->facescols+2]==1){/*Triangular base/top*/
					/*Nothing*/
				}
				else{
					_error_("not supported");
				}
			}

			/*P1 pressure*/
			vnodes = id0+3*iomodel->numberoffaces;
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(vnodes+i+1,counter+1,lid++,i,iomodel,analysis,FSpressureEnum));
				}
				counter++;
			}
			break;
		case CrouzeixRaviartEnum:
			_assert_(approximation==FSApproximationEnum);
			/*P2b velocity*/
			EdgesPartitioning(&my_edges,iomodel);
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,FSvelocityEnum));
				}
			}
			for(i=0;i<iomodel->numberofedges;i++){
				if(my_edges[i]){
					nodes->AddObject(new Node(id0+iomodel->numberofvertices+i+1,iomodel->numberofvertices+i,lid++,0,iomodel,analysis,FSvelocityEnum));
				}
			}
			id0 = id0+iomodel->numberofvertices+iomodel->numberofedges;
			if(iomodel->meshelementtype==PentaEnum){
				FacesPartitioning(&my_faces,iomodel);
				for(i=0;i<iomodel->numberoffaces;i++){
					if(iomodel->faces[i*iomodel->facescols+2]==2){/*Vertical quads*/
						if(my_faces[i]){
							node = new Node(id0+i+1,iomodel->numberofvertices+iomodel->numberofedges+i,lid++,0,iomodel,analysis,FSvelocityEnum);
							nodes->AddObject(node);
						}
					}
					else if(iomodel->faces[i*iomodel->facescols+2]==1){/*Triangular base/top*/
						/*Nothing*/
					}
					else{
						_error_("not supported");
					}
				}
				id0 = id0+iomodel->numberoffaces;
			}
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					nodes->AddObject(new Node(id0+i+1,id0-iomodel->nodecounter+i,lid++,0,iomodel,analysis,FSvelocityEnum));
				}
			}

			/*P1 DG pressure*/
			vnodes = id0+iomodel->numberofelements;
			switch(iomodel->meshelementtype){
				case TriaEnum:  elementnbv = 3; break;
				case TetraEnum: elementnbv = 4; break;
				case PentaEnum: elementnbv = 6; break;
				default:        _error_("mesh dimension not supported yet");
			}
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					for(j=0;j<elementnbv;j++){
						nodes->AddObject(new Node(vnodes+elementnbv*i+j+1,vnodes-iomodel->nodecounter+elementnbv*i+j,lid++,iomodel->elements[+elementnbv*i+j]-1,iomodel,analysis,FSpressureEnum));

					}
				}
			}
			break;
		case LACrouzeixRaviartEnum:
			_assert_(approximation==FSApproximationEnum);
			/*P2b velocity*/
			EdgesPartitioning(&my_edges,iomodel);
			for(i=0;i<iomodel->numberofvertices;i++){
				if(iomodel->my_vertices[i]){
					nodes->AddObject(new Node(id0+i+1,i,lid++,i,iomodel,analysis,FSvelocityEnum));
				}
			}
			for(i=0;i<iomodel->numberofedges;i++){
				if(my_edges[i]){
					nodes->AddObject(new Node(id0+iomodel->numberofvertices+i+1,iomodel->numberofvertices+i,lid++,0,iomodel,analysis,FSvelocityEnum));
				}
			}
			id0 = id0+iomodel->numberofvertices+iomodel->numberofedges;
			if(iomodel->meshelementtype==PentaEnum){
				FacesPartitioning(&my_faces,iomodel);
				for(i=0;i<iomodel->numberoffaces;i++){
					if(iomodel->faces[i*iomodel->facescols+2]==2){/*Vertical quads*/
						if(my_faces[i]){
							node = new Node(id0+i+1,iomodel->numberofvertices+iomodel->numberofedges+i,lid++,0,iomodel,analysis,FSvelocityEnum);
							nodes->AddObject(node);
						}
					}
					else if(iomodel->faces[i*iomodel->facescols+2]==1){/*Triangular base/top*/
						/*Nothing*/
					}
					else{
						_error_("not supported");
					}
				}
				id0 = id0+iomodel->numberoffaces;
			}
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]){
					nodes->AddObject(new Node(id0+i+1,id0-iomodel->nodecounter+i,lid++,0,iomodel,analysis,FSvelocityEnum));
				}
			}

			/*No pressure*/
			break;

		default:
			_error_("Finite element "<<EnumToStringx(finite_element)<<" not supported yet");
	}

	/*Clean up*/
	xDelete<bool>(my_faces);
	xDelete<bool>(my_edges);
	xDelete<bool>(my_nodes);
}
