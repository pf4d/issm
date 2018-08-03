/*!\file Moulin.c
 * \brief: implementation of the Moulin object
 */

/*Headers*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "shared/shared.h"
#include "../../analyses/analyses.h"
/*}}}*/

/*Element macros*/
#define NUMVERTICES   1

/*Moulin constructors and destructor*/
Moulin::Moulin(){/*{{{*/
	this->parameters=NULL;
	this->hnode=NULL;
	this->node=NULL;
	this->helement=NULL;
	this->element=NULL;
	this->hmatpar=NULL;
	this->matpar=NULL;
}
/*}}}*/
Moulin::Moulin(int id, int index, IoModel* iomodel, int in_analysis_type){ //i is the element index/*{{{*/

	int pengrid_node_id;
	int pengrid_matpar_id;
	int pengrid_element_id;

	/*Some checks if debugging activated*/
	_assert_(iomodel->singlenodetoelementconnectivity);
	_assert_(index>=0 && index<iomodel->numberofvertices);
	_assert_(id);

	/*id: */
	this->id=id;
	this->analysis_type=in_analysis_type;

	/*hooks: */
	pengrid_node_id=iomodel->nodecounter+index+1;
	pengrid_element_id=iomodel->singlenodetoelementconnectivity[index];
	_assert_(pengrid_element_id);
	pengrid_matpar_id=iomodel->numberofelements+1; //refers to the constant material parameters object

	this->hnode=new Hook(&pengrid_node_id,1);
	this->helement=new Hook(&pengrid_element_id,1);
	this->hmatpar=new Hook(&pengrid_matpar_id,1);

	//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
	this->parameters=NULL;
	this->node=NULL;
	this->element=NULL;
	this->matpar=NULL;
}
/*}}}*/
Moulin::~Moulin(){/*{{{*/
	delete hnode;
	delete helement;
	delete hmatpar;
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Moulin::copy() {/*{{{*/

	Moulin* pengrid=NULL;

	pengrid=new Moulin();

	/*copy fields: */
	pengrid->id=this->id;
	pengrid->analysis_type=this->analysis_type;

	/*point parameters: */
	pengrid->parameters=this->parameters;

	/*now deal with hooks and objects: */
	pengrid->hnode=(Hook*)this->hnode->copy();
	pengrid->hmatpar=(Hook*)this->hmatpar->copy();
	pengrid->helement=(Hook*)this->helement->copy();

	/*corresponding fields*/
	pengrid->node  =(Node*)pengrid->hnode->delivers();
	pengrid->matpar =(Matpar*)pengrid->hmatpar->delivers();
	pengrid->element=(Element*)pengrid->helement->delivers();

	return pengrid;
}
/*}}}*/
void    Moulin::DeepEcho(void){/*{{{*/

	_printf_("Moulin:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	hnode->DeepEcho();
	helement->DeepEcho();
	hmatpar->DeepEcho();
	_printf_("   parameters\n");
	parameters->DeepEcho();
}
/*}}}*/
void    Moulin::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int     Moulin::Id(void){ return id; }/*{{{*/
/*}}}*/
void    Moulin::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	_assert_(this);

	/*ok, marshall operations: */
	MARSHALLING_ENUM(MoulinEnum);
	MARSHALLING(id);
	MARSHALLING(analysis_type);

	if(marshall_direction==MARSHALLING_BACKWARD){
		this->hnode      = new Hook();
		this->helement   = new Hook();
		this->hmatpar    = new Hook();
	}

	this->hnode->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->helement->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->hmatpar->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);

	/*corresponding fields*/
	node   =(Node*)this->hnode->delivers();
	matpar =(Matpar*)this->hmatpar->delivers();
	element=(Element*)this->helement->delivers();
}
/*}}}*/
int     Moulin::ObjectEnum(void){/*{{{*/

	return MoulinEnum;
}
/*}}}*/

/*Load virtual functions definitions:*/
void  Moulin::Configure(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

	/*Take care of hooking up all objects for this load, ie links the objects in the hooks to their respective 
	 * datasets, using internal ids and offsets hidden in hooks: */
	hnode->configure(nodesin);
	helement->configure(elementsin);
	hmatpar->configure(materialsin);

	/*Get corresponding fields*/
	node=(Node*)hnode->delivers();
	element=(Element*)helement->delivers();
	matpar=(Matpar*)hmatpar->delivers();

	/*point parameters to real dataset: */
	this->parameters=parametersin;
}
/*}}}*/
void  Moulin::CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/

	/*No loads applied, do nothing: */
	return;

}
/*}}}*/
void  Moulin::CreatePVector(Vector<IssmDouble>* pf){/*{{{*/

	ElementVector* pe=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		
	case HydrologySommersAnalysisEnum:
		pe = this->CreatePVectorHydrologySommers();
		break;
	case HydrologyDCInefficientAnalysisEnum:
		pe = CreatePVectorHydrologyDCInefficient();
		break;
	case HydrologyDCEfficientAnalysisEnum:
		pe = CreatePVectorHydrologyDCEfficient();
		break;
	default:
		_error_("Don't know why we should be here");
		/*No loads applied, do nothing: */
		return;
	}
	if(pe){
		pe->AddToGlobal(pf);
		delete pe;
	}

}
/*}}}*/
void  Moulin::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(node);

	lidlist[0]=node->Lid();
}
/*}}}*/
void  Moulin::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(node);

	sidlist[0]=node->Sid();
}
/*}}}*/
int   Moulin::GetNumberOfNodes(void){/*{{{*/

	return NUMVERTICES;
}
/*}}}*/
bool  Moulin::InAnalysis(int in_analysis_type){/*{{{*/
	if (in_analysis_type==this->analysis_type)return true;
	else return false;
}
/*}}}*/
bool  Moulin::IsPenalty(void){/*{{{*/
	return true;
}
/*}}}*/
void  Moulin::PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs,IssmDouble kmax){/*{{{*/

	/*Don't do anything for now*/

}
/*}}}*/
void  Moulin::PenaltyCreatePVector(Vector<IssmDouble>* pf,IssmDouble kmax){/*{{{*/

	/*Don't do anything for now*/
}
/*}}}*/
void  Moulin::ResetHooks(){/*{{{*/

	this->node=NULL;
	this->element=NULL;
	this->matpar=NULL;
	this->parameters=NULL;

	/*Get Element type*/
	this->hnode->reset();
	this->helement->reset();
	this->hmatpar->reset();

}
/*}}}*/
void  Moulin::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Moulin::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int set1_enum,int set2_enum){/*{{{*/

	/*Output */
	int d_nz = 0;
	int o_nz = 0;

	if(!flags[this->node->Lid()]){

		/*flag current node so that no other element processes it*/
		flags[this->node->Lid()]=true;

		int counter=0;
		while(flagsindices[counter]>=0) counter++;
		flagsindices[counter]=this->node->Lid();

		/*if node is clone, we have an off-diagonal non-zero, else it is a diagonal non-zero*/
		switch(set2_enum){
			case FsetEnum:
				if(node->indexing.fsize){
					if(this->node->IsClone())
					 o_nz += 1;
					else
					 d_nz += 1;
				}
				break;
			case GsetEnum:
				if(node->indexing.gsize){
					if(this->node->IsClone())
					 o_nz += 1;
					else
					 d_nz += 1;
				}
				break;
			case SsetEnum:
				if(node->indexing.ssize){
					if(this->node->IsClone())
					 o_nz += 1;
					else
					 d_nz += 1;
				}
				break;
			default: _error_("not supported");
		}
	}

	/*Assign output pointers: */
	*pd_nz=d_nz;
	*po_nz=o_nz;
}
/*}}}*/

/*Update virtual functions definitions:*/
void  Moulin::InputUpdateFromConstant(IssmDouble constant, int name){/*{{{*/
	/*Nothing*/
}
/*}}}*/
void  Moulin::InputUpdateFromConstant(int constant, int name){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Moulin::InputUpdateFromConstant(bool constant, int name){/*{{{*/

	/*Don't do anything for now*/
}
/*}}}*/
void  Moulin::InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrows, int ncols, int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Moulin::InputUpdateFromVector(IssmDouble* vector, int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Moulin::InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/

ElementVector* Moulin::CreatePVectorHydrologySommers(void){/*{{{*/

	/*If this node is not the master node (belongs to another partition of the
	 * mesh), don't add the moulin input a second time*/
	if(node->IsClone()) return NULL;

	IssmDouble moulin_load;

	/*Initialize Element matrix*/
	ElementVector* pe=new ElementVector(&node,1,this->parameters);

	this->element->GetInputValue(&moulin_load,node,HydrologyMoulinInputEnum);
	pe->values[0]=moulin_load;

	/*Clean up and return*/
	return pe;
}
/*}}}*/
ElementVector* Moulin::CreatePVectorHydrologyDCInefficient(void){/*{{{*/

	/*If this node is not the master node (belongs to another partition of the
	 * mesh), don't add the moulin input a second time*/
	if(node->IsClone()) return NULL;
	bool isefficientlayer;
	IssmDouble moulin_load,dt;
	IssmDouble epl_active;

	/*Initialize Element matrix*/
	ElementVector* pe=new ElementVector(&node,1,this->parameters);

	this->element->GetInputValue(&moulin_load,node,HydrologydcBasalMoulinInputEnum);
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	parameters->FindParam(&isefficientlayer,HydrologydcIsefficientlayerEnum);
	// Test version input in EPL when active
	if(isefficientlayer){
		this->element->GetInputValue(&epl_active,node,HydrologydcMaskEplactiveNodeEnum);
		if(reCast<bool>(epl_active)){
			pe->values[0]=moulin_load*0.0;
		}
		else{
			pe->values[0]=moulin_load*dt;
		}
	}
	else{
		pe->values[0]=moulin_load*dt;
	}
	/*Clean up and return*/
	return pe;
 }
/*}}}*/
ElementVector* Moulin::CreatePVectorHydrologyDCEfficient(void){/*{{{*/

	/*If this node is not the master node (belongs to another partition of the
	 * mesh), don't add the moulin input a second time*/
	if(node->IsClone()) return NULL;
	if(!this->node->IsActive()) return NULL;
	IssmDouble moulin_load,dt;
	ElementVector* pe=new ElementVector(&node,1,this->parameters);
	
	this->element->GetInputValue(&moulin_load,node,HydrologydcBasalMoulinInputEnum);
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	
	pe->values[0]=moulin_load*dt;
	/*Clean up and return*/
	return pe;
}
/*}}}*/
