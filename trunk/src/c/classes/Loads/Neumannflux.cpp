/*!\file Neumannflux.c
 * \brief: implementation of the Neumannflux object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "shared/shared.h"
#include "../classes.h"
/*}}}*/	

/*Load macros*/
#define NUMVERTICES 2
#define NUMNODES_BOUNDARY 2

/*Neumannflux constructors and destructor*/
Neumannflux::Neumannflux(){/*{{{*/
	this->parameters = NULL;
	this->helement   = NULL;
	this->element    = NULL;
	this->hnodes     = NULL;
	this->hvertices  = NULL;
	this->nodes      = NULL;
}
/*}}}*/
Neumannflux::Neumannflux(int neumannflux_id,int i,IoModel* iomodel,int* segments,int in_analysis_type){/*{{{*/


	/*Some sanity checks*/
	_assert_(segments);

	/*neumannflux constructor data: */
	int neumannflux_elem_id;
	int neumannflux_vertex_ids[2];
	int neumannflux_node_ids[2];

	/*1: Get vertices ids*/
	neumannflux_vertex_ids[0]=segments[3*i+0];
	neumannflux_vertex_ids[1]=segments[3*i+1];

	/*2: Get the ids of the nodes*/
	neumannflux_node_ids[0]=iomodel->nodecounter+neumannflux_vertex_ids[0];
	neumannflux_node_ids[1]=iomodel->nodecounter+neumannflux_vertex_ids[1];

	/*Get element id*/
	neumannflux_elem_id = segments[3*i+2];

	/*Ok, we have everything to build the object: */
	this->id=neumannflux_id;
	this->analysis_type=in_analysis_type;

	/*Hooks: */
	this->hnodes    =new Hook(&neumannflux_node_ids[0],2);
	this->hvertices =new Hook(&neumannflux_vertex_ids[0],2);
	this->helement  =new Hook(&neumannflux_elem_id,1);

	//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
	this->parameters=NULL;
	this->element=NULL;
	this->nodes=NULL;
}
/*}}}*/
Neumannflux::~Neumannflux(){/*{{{*/
	this->parameters=NULL;
	delete helement;
	delete hnodes;
	delete hvertices;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Neumannflux::copy() {/*{{{*/

	Neumannflux* neumannflux=NULL;

	neumannflux=new Neumannflux();

	/*copy fields: */
	neumannflux->id=this->id;
	neumannflux->analysis_type=this->analysis_type;

	/*point parameters: */
	neumannflux->parameters=this->parameters;

	/*now deal with hooks and objects: */
	neumannflux->hnodes    = (Hook*)this->hnodes->copy();
	neumannflux->hvertices = (Hook*)this->hvertices->copy();
	neumannflux->helement  = (Hook*)this->helement->copy();

	/*corresponding fields*/
	neumannflux->nodes    = (Node**)neumannflux->hnodes->deliverp();
	neumannflux->vertices = (Vertex**)neumannflux->hvertices->deliverp();
	neumannflux->element  = (Element*)neumannflux->helement->delivers();

	return neumannflux;
}
/*}}}*/
void    Neumannflux::DeepEcho(void){/*{{{*/

	_printf_("Neumannflux:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	hnodes->DeepEcho();
	hvertices->DeepEcho();
	helement->DeepEcho();
	_printf_("   parameters\n");
	if(parameters)
	 parameters->DeepEcho();
	else
	 _printf_("      NULL\n");
}		
/*}}}*/
void    Neumannflux::Echo(void){/*{{{*/
	_printf_("Neumannflux:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	hnodes->Echo();
	hvertices->Echo();
	helement->Echo();
	_printf_("   parameters: " << parameters << "\n");
}
/*}}}*/
int     Neumannflux::Id(void){/*{{{*/
	return id;
}
/*}}}*/
void    Neumannflux::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	_assert_(this);

	/*ok, marshall operations: */
	MARSHALLING_ENUM(NeumannfluxEnum);
	MARSHALLING(id);
	MARSHALLING(analysis_type);

	if(marshall_direction==MARSHALLING_BACKWARD){
		this->hnodes      = new Hook();
		this->hvertices   = new Hook();
		this->helement    = new Hook();
	}

	this->hnodes->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->helement->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->hvertices->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);

	/*corresponding fields*/
	nodes    =(Node**)this->hnodes->deliverp();
	vertices =(Vertex**)this->hvertices->deliverp();
	element  =(Element*)this->helement->delivers();

}
/*}}}*/
int     Neumannflux::ObjectEnum(void){/*{{{*/

	return NeumannfluxEnum;

}
/*}}}*/

/*Load virtual functions definitions:*/
void  Neumannflux::Configure(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective 
	 * datasets, using internal ids and offsets hidden in hooks: */
	hnodes->configure((DataSet*)nodesin);
	hvertices->configure((DataSet*)verticesin);
	helement->configure((DataSet*)elementsin);

	/*Initialize hooked fields*/
	this->nodes    = (Node**)hnodes->deliverp();
	this->vertices = (Vertex**)hvertices->deliverp();
	this->element  = (Element*)helement->delivers();

	/*point parameters to real dataset: */
	this->parameters=parametersin;
}
/*}}}*/
void  Neumannflux::CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/

	/*recover some parameters*/
	ElementMatrix* Ke=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	/*Just branch to the correct element stiffness matrix generator, according to the type of analysis we are carrying out: */
	switch(analysis_type){
		case HydrologySommersAnalysisEnum:
			/*Nothing!*/
			break;
		default:
			_error_("analysis " << analysis_type << " (" << EnumToStringx(analysis_type) << ") not supported yet");
	}

	/*Add to global matrix*/
	if(Ke){
		Ke->AddToGlobal(Kff,Kfs);
		delete Ke;
	}

}
/*}}}*/
void  Neumannflux::CreatePVector(Vector<IssmDouble>* pf){/*{{{*/

	/*recover some parameters*/
	ElementVector* pe=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case HydrologySommersAnalysisEnum:
			pe=CreatePVectorHydrologySommers();
			break;
		default:
			_error_("analysis " << analysis_type << " (" << EnumToStringx(analysis_type) << ") not supported yet");
	}

	/*Add to global matrix*/
	if(pe){
		pe->AddToGlobal(pf);
		delete pe;
	}

}
/*}}}*/
void  Neumannflux::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(nodes);

	for(int i=0;i<NUMNODES_BOUNDARY;i++) lidlist[i]=nodes[i]->Lid();
}
/*}}}*/
void  Neumannflux::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(nodes);

	for(int i=0;i<NUMNODES_BOUNDARY;i++) sidlist[i]=nodes[i]->Sid();
}
/*}}}*/
int   Neumannflux::GetNumberOfNodes(void){/*{{{*/

	return NUMNODES_BOUNDARY;
}
/*}}}*/
bool  Neumannflux::InAnalysis(int in_analysis_type){/*{{{*/
	if (in_analysis_type==this->analysis_type) return true;
	else return false;
}
/*}}}*/
bool  Neumannflux::IsPenalty(void){/*{{{*/
	return false;
}
/*}}}*/
void  Neumannflux::PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs,IssmDouble kmax){/*{{{*/

	/*No stiffness loads applied, do nothing: */
	return;

}
/*}}}*/
void  Neumannflux::PenaltyCreatePVector(Vector<IssmDouble>* pf,IssmDouble kmax){/*{{{*/

	/*No penalty loads applied, do nothing: */
	return;

}
/*}}}*/
void  Neumannflux::ResetHooks(){/*{{{*/

	this->nodes=NULL;
	this->vertices=NULL;
	this->element=NULL;
	this->parameters=NULL;

	/*Get Element type*/
	this->hnodes->reset();
	this->hvertices->reset();
	this->helement->reset();

}
/*}}}*/
void  Neumannflux::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Neumannflux::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int set1_enum,int set2_enum){/*{{{*/

	/*Output */
	int d_nz = 0;
	int o_nz = 0;

	/*Loop over all nodes*/
	for(int i=0;i<this->GetNumberOfNodes();i++){

		if(!flags[this->nodes[i]->Lid()]){

			/*flag current node so that no other element processes it*/
			flags[this->nodes[i]->Lid()]=true;

			int counter=0;
			while(flagsindices[counter]>=0) counter++;
			flagsindices[counter]=this->nodes[i]->Lid();

			/*if node is clone, we have an off-diagonal non-zero, else it is a diagonal non-zero*/
			switch(set2_enum){
				case FsetEnum:
					if(nodes[i]->indexing.fsize){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				case GsetEnum:
					if(nodes[i]->indexing.gsize){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				case SsetEnum:
					if(nodes[i]->indexing.ssize){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				default: _error_("not supported");
			}
		}
	}

	/*Assign output pointers: */
	*pd_nz=d_nz;
	*po_nz=o_nz;
}
/*}}}*/

/*Neumannflux management*/
ElementVector* Neumannflux::CreatePVectorHydrologySommers(void){/*{{{*/

	/* constants*/
	const int numdof=2;

	/* Intermediaries*/
	IssmDouble Jdet,flux;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble basis[numdof];

	/*Initialize Load Vector and return if necessary*/
	Tria*  tria=(Tria*)element;
	_assert_(tria->FiniteElement()==P1Enum); 
	if(!tria->IsIceInElement() || tria->IsFloating()) return NULL;

	/*Initialize Element vector and other vectors*/
	ElementVector* pe=new ElementVector(nodes,NUMNODES_BOUNDARY,this->parameters);

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* flux_input = tria->inputs->GetInput(HydrologyNeumannfluxEnum);  _assert_(flux_input); 

	/*Check wether it is an inflow or outflow BC (0 is the middle of the segment)*/
	int index1=tria->GetNodeIndex(nodes[0]);
	int index2=tria->GetNodeIndex(nodes[1]);

	/* Start  looping on the number of gaussian points: */
	GaussTria* gauss=new GaussTria(index1,index2,2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		tria->GetSegmentNodalFunctions(&basis[0],gauss,index1,index2,tria->FiniteElement());
		flux_input->GetInputValue(&flux,gauss);

		for(int i=0;i<numdof;i++) pe->values[i] += gauss->weight*Jdet*flux*basis[i];
	}

	/*Clean up and return*/
	delete gauss;
	return pe;
}
/*}}}*/
