/*!\file Riftfront.cpp
 * \brief: implementation of the Riftfront object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "shared/shared.h"
#include "modules/ModelProcessorx/ModelProcessorx.h"
#include "../classes.h"
/*}}}*/

/*Element macros*/
#define NUMVERTICES 2

/*Riftfront constructors and destructor*/
Riftfront::Riftfront(){/*{{{*/
	this->parameters=NULL;
	this->hnodes=NULL;
	this->helements=NULL;
	this->hmatpar=NULL;
	this->nodes=NULL;
	this->elements=NULL;
	this->matpar=NULL;
}
/*}}}*/
Riftfront::Riftfront(int riftfront_id,int i, IoModel* iomodel,int riftfront_analysis_type){/*{{{*/

	/*data: */
	const int RIFTINFOSIZE = 12;
	int    riftfront_node_ids[2];
	int    riftfront_elem_ids[2];
	int    riftfront_matpar_id;
	IssmDouble riftfront_friction;
	IssmDouble riftfront_fractionincrement;
	int    penalty_lock;

	/*intermediary: */
	int el1    ,el2;
	int node1  ,node2;

	/*Fetch parameters: */
	iomodel->FindConstant(&penalty_lock,"md.stressbalance.rift_penalty_lock");

	/*Ok, retrieve all the data needed to add a penalty between the two nodes: */
	el1=reCast<int,IssmDouble>(*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+2));
	el2=reCast<int,IssmDouble>(*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+3)) ;

	node1=reCast<int,IssmDouble>(*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+0));
	node2=reCast<int,IssmDouble>(*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+1));

	/*id: */
	this->id=riftfront_id;
	this->analysis_type=riftfront_analysis_type;

	/*hooks: */
	riftfront_node_ids[0]=iomodel->nodecounter+node1;
	riftfront_node_ids[1]=iomodel->nodecounter+node2;
	riftfront_elem_ids[0]=el1;
	riftfront_elem_ids[1]=el2;
	riftfront_matpar_id=iomodel->numberofelements+1; //matlab indexing

	/*Hooks: */
	this->hnodes=new Hook(riftfront_node_ids,2);
	this->helements=new Hook(riftfront_elem_ids,2);
	this->hmatpar=new Hook(&riftfront_matpar_id,1);

	/*computational parameters: */
	this->active=0;
	this->frozen=0;
	this->counter=0;
	this->prestable=0;
	this->penalty_lock=penalty_lock;
	this->material_converged=0;
	this->normal[0]=*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+4);
	this->normal[1]=*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+5);
	this->length=*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+6);
	this->fraction=*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+9);
	this->state=reCast<int,IssmDouble>(*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+11));

	//intialize properties
	this->type=SegmentRiftfrontEnum;
	this->fill = IoRiftfillToEnum(reCast<int,IssmDouble>(*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+7)));
	this->friction=*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+8);
	this->fractionincrement=*(iomodel->Data("md.rifts.riftstruct")+RIFTINFOSIZE*i+10);
	this->shelf=reCast<bool,IssmDouble>(iomodel->Data("md.mask.groundedice_levelset")[node1-1]<0.);

	//parameters and hooked fields: we still can't point to them, they may not even exist. Configure will handle this.
	this->parameters=NULL;
	this->nodes= NULL;
	this->elements= NULL;
	this->matpar= NULL;

}
/*}}}*/
Riftfront::~Riftfront(){/*{{{*/
	this->parameters=NULL;
	delete hnodes;
	delete helements;
	delete hmatpar;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Riftfront::copy() {/*{{{*/

	Riftfront* riftfront=NULL;

	riftfront=new Riftfront();

	/*copy fields: */
	riftfront->id=this->id;
	riftfront->analysis_type=this->analysis_type;
	riftfront->type=this->type;
	riftfront->fill=this->fill;
	riftfront->friction=this->friction;
	riftfront->fractionincrement=this->fractionincrement;
	riftfront->shelf=this->shelf;

	/*point parameters: */
	riftfront->parameters=this->parameters;

	/*now deal with hooks and objects: */
	riftfront->hnodes=(Hook*)this->hnodes->copy();
	riftfront->helements=(Hook*)this->helements->copy();
	riftfront->hmatpar=(Hook*)this->hmatpar->copy();

	/*corresponding fields*/
	riftfront->nodes   =(Node**)riftfront->hnodes->deliverp();
	riftfront->elements=(Element**)riftfront->helements->deliverp();
	riftfront->matpar  =(Matpar*)riftfront->hmatpar->delivers();

	/*internal data: */
	riftfront->penalty_lock=this->penalty_lock;
	riftfront->active=this->active;
	riftfront->frozen=this->frozen;
	riftfront->state=this->state;
	riftfront->counter=this->counter;
	riftfront->prestable=this->prestable;
	riftfront->material_converged=this->material_converged;
	riftfront->normal[0]=this->normal[0];
	riftfront->normal[1]=this->normal[1];
	riftfront->length=this->length;
	riftfront->fraction=this->fraction;

	return riftfront;

}
/*}}}*/
void    Riftfront::DeepEcho(void){/*{{{*/

	_printf_("Riftfront:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	hnodes->DeepEcho();
	helements->DeepEcho();
	hmatpar->DeepEcho();
	_printf_("   parameters\n");
	if(parameters)parameters->DeepEcho();
}
/*}}}*/
void    Riftfront::Echo(void){/*{{{*/

	_printf_("Riftfront:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	_printf_("   hnodes: " << hnodes << "\n");
	_printf_("   helements: " << helements << "\n");
	_printf_("   hmatpar: " << hmatpar << "\n");
	_printf_("   parameters: " << parameters << "\n");
	_printf_("   internal parameters: \n");
	_printf_("   normal: " << normal[0] << "|" << normal[1] << "\n");
	_printf_("   length: " << length << "\n");
	_printf_("   penalty_lock: " << penalty_lock << "\n");
	_printf_("   active: " <<(active ? "true":"false") << "\n");
	_printf_("   counter: " << counter << "\n");
	_printf_("   prestable: " << (prestable ? "true":"false") << "\n");
	_printf_("   material_converged: " << (material_converged ? "true":"false") << "\n");
	_printf_("   fill: " << fill << "\n");
	_printf_("   friction: " << friction << "\n");
	_printf_("   fraction: " << fraction << "\n");
	_printf_("   fractionincrement: " << fractionincrement << "\n");
	_printf_("   state: " << state << "\n");
	_printf_("   frozen: " << (frozen ? "true":"false") << "\n");

}
/*}}}*/
int     Riftfront::Id(void){ return id; }/*{{{*/
/*}}}*/
void    Riftfront::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	_assert_(this);

	/*ok, marshall operations: */
	MARSHALLING_ENUM(RiftfrontEnum);
	MARSHALLING(id);
	MARSHALLING(analysis_type);
	MARSHALLING(type);
	MARSHALLING(fill);
	MARSHALLING(friction);
	MARSHALLING(fractionincrement);
	MARSHALLING(shelf);

	if(marshall_direction==MARSHALLING_BACKWARD){
		this->hnodes      = new Hook();
		this->hmatpar     = new Hook();
		this->helements   = new Hook();
	}

	this->hnodes->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->hmatpar->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->helements->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);

	/*corresponding fields*/
	nodes     =(Node**)this->hnodes->deliverp();
	matpar    =(Matpar*)this->hmatpar->delivers();
	elements  =(Element**)this->helements->deliverp();

	MARSHALLING(penalty_lock);
	MARSHALLING(active);
	MARSHALLING(frozen);
	MARSHALLING(state);
	MARSHALLING(counter);
	MARSHALLING(prestable);
	MARSHALLING(material_converged);
	MARSHALLING(normal[0]);
	MARSHALLING(normal[1]);
	MARSHALLING(length);
	MARSHALLING(fraction);

}
/*}}}*/
int     Riftfront::ObjectEnum(void){/*{{{*/

	return RiftfrontEnum;

}
/*}}}*/

/*Update virtual functions definitions:*/
void  Riftfront::InputUpdateFromConstant(bool constant,int name){/*{{{*/
}
/*}}}*/
void  Riftfront::InputUpdateFromConstant(IssmDouble constant,int name){/*{{{*/

}
/*}}}*/
void  Riftfront::InputUpdateFromVector(IssmDouble* vector, int name, int type){/*{{{*/

	/*Nothing to update*/

}
/*}}}*/

/*Load virtual functions definitions:*/
void  Riftfront::Configure(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/	

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective 
	 * datasets, using internal ids and offsets hidden in hooks: */
	hnodes->configure(nodesin);
	helements->configure(elementsin);
	hmatpar->configure(materialsin);

	/*Initialize hooked fields*/
	this->nodes   =(Node**)hnodes->deliverp();
	this->elements=(Element**)helements->deliverp();
	this->matpar  =(Matpar*)hmatpar->delivers();

	/*point parameters to real dataset: */
	this->parameters=parametersin;

}
/*}}}*/
void  Riftfront::CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/
	/*do nothing: */
	return;
}
/*}}}*/
void  Riftfront::CreatePVector(Vector<IssmDouble>* pf){/*{{{*/
	/*do nothing: */
	return;
}
/*}}}*/
void  Riftfront::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(nodes);

	for(int i=0;i<NUMVERTICES;i++) lidlist[i]=nodes[i]->Lid();
}
/*}}}*/
void  Riftfront::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(nodes);

	for(int i=0;i<NUMVERTICES;i++) sidlist[i]=nodes[i]->Sid();
}
/*}}}*/
int   Riftfront::GetNumberOfNodes(void){/*{{{*/

	return NUMVERTICES;
}
/*}}}*/
bool  Riftfront::InAnalysis(int in_analysis_type){/*{{{*/
	if (in_analysis_type==this->analysis_type) return true;
	else return false;
}
/*}}}*/
bool  Riftfront::IsPenalty(void){/*{{{*/
	return true;
}
/*}}}*/
void  Riftfront::PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs,IssmDouble kmax){/*{{{*/

	/*Retrieve parameters: */
	ElementMatrix* Ke=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case StressbalanceAnalysisEnum:
			Ke=PenaltyCreateKMatrixStressbalanceHoriz(kmax);
			break;
		case AdjointHorizAnalysisEnum:
			Ke=PenaltyCreateKMatrixStressbalanceHoriz(kmax);
			break;
		default:
			_error_("analysis " << analysis_type << " (" << EnumToStringx(analysis_type) << ") not supported yet");
	}

	/*Add to global Vector*/
	if(Ke){
		Ke->AddToGlobal(Kff,Kfs);
		delete Ke;
	}
}
/*}}}*/
void  Riftfront::PenaltyCreatePVector(Vector<IssmDouble>* pf,IssmDouble kmax){/*{{{*/

	/*Retrieve parameters: */
	ElementVector* pe=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case StressbalanceAnalysisEnum:
			pe=PenaltyCreatePVectorStressbalanceHoriz(kmax);
			break;
		case AdjointHorizAnalysisEnum:
			/*No penalty applied on load vector*/
			break;
		default:
			_error_("analysis " << analysis_type << " (" << EnumToStringx(analysis_type) << ") not supported yet");
	}

	/*Add to global Vector*/
	if(pe){
		pe->AddToGlobal(pf);
		delete pe;
	}
}
/*}}}*/
void  Riftfront::ResetHooks(){/*{{{*/

	this->nodes=NULL;
	this->elements=NULL;
	this->matpar=NULL;
	this->parameters=NULL;

	/*Get Element type*/
	this->hnodes->reset();
	this->helements->reset();
	this->hmatpar->reset();

}
/*}}}*/
void  Riftfront::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Riftfront::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int set1_enum,int set2_enum){/*{{{*/

	/*Output */
	int d_nz = 0;
	int o_nz = 0;

	/*Loop over all nodes*/
	for(int i=0;i<NUMVERTICES;i++){

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

/*Riftfront numerics*/
ElementMatrix* Riftfront::PenaltyCreateKMatrixStressbalanceHoriz(IssmDouble kmax){/*{{{*/

	const int   numdof = NDOF2*NUMVERTICES;
	IssmDouble  thickness;
	IssmDouble  h[2];
	IssmDouble  penalty_offset;

	/*Objects: */
	Tria       *tria1               = NULL;
	Tria       *tria2               = NULL;

	/*enum of element? */
	if(elements[0]->ObjectEnum()!=TriaEnum)_error_("only Tria element allowed for Riftfront load!");
	tria1=(Tria*)elements[0];
	tria2=(Tria*)elements[1];

	/*Initialize Element Matrix*/
	if(!this->active) return NULL;
	ElementMatrix* Ke=new ElementMatrix(nodes,NUMVERTICES,this->parameters);

	/*Get some parameters: */
	this->parameters->FindParam(&penalty_offset,StressbalancePenaltyFactorEnum);
	tria1->GetInputValue(&h[0],nodes[0],ThicknessEnum);
	tria2->GetInputValue(&h[1],nodes[1],ThicknessEnum);
	if (h[0]!=h[1])_error_("different thicknesses not supported for rift fronts");
	thickness=h[0];

	/*There is contact, we need to constrain the normal velocities (zero penetration), and the 
	 *contact slip friction. */

	/*From Peter Wriggers book (Computational Contact Mechanics, p191): */
	Ke->values[0*numdof+0]+= +pow(normal[0],2)*kmax*pow(10,penalty_offset);
	Ke->values[0*numdof+1]+= +normal[0]*normal[1]*kmax*pow(10,penalty_offset);
	Ke->values[0*numdof+2]+= -pow(normal[0],2)*kmax*pow(10,penalty_offset);
	Ke->values[0*numdof+3]+= -normal[0]*normal[1]*kmax*pow(10,penalty_offset);

	Ke->values[1*numdof+0]+= +normal[0]*normal[1]*kmax*pow(10,penalty_offset);
	Ke->values[1*numdof+1]+= +pow(normal[1],2)*kmax*pow(10,penalty_offset);
	Ke->values[1*numdof+2]+= -normal[0]*normal[1]*kmax*pow(10,penalty_offset);
	Ke->values[1*numdof+3]+= -pow(normal[1],2)*kmax*pow(10,penalty_offset);

	Ke->values[2*numdof+0]+= -pow(normal[0],2)*kmax*pow(10,penalty_offset);
	Ke->values[2*numdof+1]+= -normal[0]*normal[1]*kmax*pow(10,penalty_offset);
	Ke->values[2*numdof+2]+= +pow(normal[0],2)*kmax*pow(10,penalty_offset);
	Ke->values[2*numdof+3]+= +normal[0]*normal[1]*kmax*pow(10,penalty_offset);

	Ke->values[3*numdof+0]+= -normal[0]*normal[1]*kmax*pow(10,penalty_offset);
	Ke->values[3*numdof+1]+= -pow(normal[1],2)*kmax*pow(10,penalty_offset);
	Ke->values[3*numdof+2]+= +normal[0]*normal[1]*kmax*pow(10,penalty_offset);
	Ke->values[3*numdof+3]+= +pow(normal[1],2)*kmax*pow(10,penalty_offset);

	/*Now take care of the friction: of type sigma=frictiontangent_velocity2-tangent_velocity1)*/

	Ke->values[0*numdof+0]+= +pow(normal[1],2)*thickness*length*friction;
	Ke->values[0*numdof+1]+= -normal[0]*normal[1]*thickness*length*friction;
	Ke->values[0*numdof+2]+= -pow(normal[1],2)*thickness*length*friction;
	Ke->values[0*numdof+3]+= +normal[0]*normal[1]*thickness*length*friction;

	Ke->values[1*numdof+0]+= -normal[0]*normal[1]*thickness*length*friction;
	Ke->values[1*numdof+1]+= +pow(normal[0],2)*thickness*length*friction;
	Ke->values[1*numdof+2]+= +normal[0]*normal[1]*thickness*length*friction;
	Ke->values[1*numdof+3]+= -pow(normal[0],2)*thickness*length*friction;

	Ke->values[2*numdof+0]+= -pow(normal[1],2)*thickness*length*friction;
	Ke->values[2*numdof+1]+= +normal[0]*normal[1]*thickness*length*friction;
	Ke->values[2*numdof+2]+= +pow(normal[1],2)*thickness*length*friction;
	Ke->values[2*numdof+3]+= -normal[0]*normal[1]*thickness*length*friction;

	Ke->values[3*numdof+0]+= +normal[0]*normal[1]*thickness*length*friction;
	Ke->values[3*numdof+1]+= -pow(normal[0],2)*thickness*length*friction;
	Ke->values[3*numdof+2]+= -normal[0]*normal[1]*thickness*length*friction;
	Ke->values[3*numdof+3]+= +pow(normal[0],2)*thickness*length*friction;

	/*Clean up and return*/
	return Ke;
}
/*}}}*/
ElementVector* Riftfront::PenaltyCreatePVectorStressbalanceHoriz(IssmDouble kmax){/*{{{*/

	int        j;
	IssmDouble rho_ice;
	IssmDouble rho_water;
	IssmDouble gravity;
	IssmDouble thickness;
	IssmDouble h[2];
	IssmDouble bed;
	IssmDouble b[2];
	IssmDouble pressure;
	IssmDouble pressure_litho;
	IssmDouble pressure_air;
	IssmDouble pressure_melange;
	IssmDouble pressure_water;

	/*Objects: */
	Tria *tria1 = NULL;
	Tria *tria2 = NULL;

	/*enum of element? */
	if(elements[0]->ObjectEnum()!=TriaEnum)_error_("only Tria element allowed for Riftfront load!");
	tria1=(Tria*)elements[0];
	tria2=(Tria*)elements[1];

	/*Initialize Element Matrix*/
	if(this->active) return NULL; /*The penalty is active. No loads implied here.*/
	ElementVector* pe=new ElementVector(nodes,NUMVERTICES,this->parameters);

	/*Get some inputs: */
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	rho_water=matpar->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	gravity=matpar->GetMaterialParameter(ConstantsGEnum);
	tria1->GetInputValue(&h[0],nodes[0],ThicknessEnum);
	tria2->GetInputValue(&h[1],nodes[1],ThicknessEnum);
	if (h[0]!=h[1])_error_("different thicknesses not supported for rift fronts");
	thickness=h[0];
	tria1->GetInputValue(&b[0],nodes[0],BaseEnum);
	tria2->GetInputValue(&b[1],nodes[1],BaseEnum);
	if (b[0]!=b[1])_error_("different beds not supported for rift fronts");
	bed=b[0];

	/*Ok, this rift is opening. We should put loads on both sides of the rift flanks. Because we are dealing with contact mechanics, 
	 * and we want to avoid zigzagging of the loads, we want lump the loads onto nodes, not onto surfaces between nodes.:*/

	/*Ok, to compute the pressure, we are going to need material properties, thickness and bed for the two nodes. We assume those properties to 
	 * be the same across the rift.: */

	/*Ok, now compute the pressure (in norm) that is being applied to the flanks, depending on the type of fill: */
	if(fill==WaterEnum){
		if(shelf){
			/*We are on an ice shelf, hydrostatic equilibrium is used to determine the pressure for water fill: */
			pressure=rho_ice*gravity*pow(thickness,2)/2.- rho_water*gravity*pow(bed,2)/2.; 
		}
		else{
			//We are on an icesheet, we assume the water column fills the entire front: */
			pressure=rho_ice*gravity*pow(thickness,2)/2.- rho_water*gravity*pow(thickness,2)/2.; 
		}
	}
	else if(fill==AirEnum){
		pressure=rho_ice*gravity*pow(thickness,2)/2.;   //icefront on an ice sheet, pressure imbalance ice vs air.
	}
	else if(fill==IceEnum){ //icefront finding itself against another icefront (pressure imbalance is fully compensated, ice vs ice)
		pressure=0;
	}
	else if(fill==MelangeEnum){ //icefront finding itself against another icefront (pressure imbalance is fully compensated, ice vs ice)

		if(!shelf) _error_("fill type " << fill << " not supported on ice sheets yet.");

		pressure_litho=rho_ice*gravity*pow(thickness,2)/2.;
		pressure_air=0;
		pressure_melange=rho_ice*gravity*pow(fraction*thickness,2)/2.;
		pressure_water=1.0/2.0*rho_water*gravity*  ( pow(bed,2.0)-pow(rho_ice/rho_water*fraction*thickness,2.0) );

		pressure=pressure_litho-pressure_air-pressure_melange-pressure_water;
	}
	else{
		_error_("fill type " << fill << " not supported yet.");
	}

	/*Ok, add contribution to first node, along the normal i==0: */
	for (j=0;j<2;j++){
		pe->values[j]+=pressure*normal[j]*length;
	}

	/*Add contribution to second node, along the opposite normal: i==1 */
	for (j=0;j<2;j++){
		pe->values[2+j]+= -pressure*normal[j]*length;
	}	

	/*Clean up and return*/
	return pe;
}
/*}}}*/
#define _ZIGZAGCOUNTER_
int    Riftfront::Constrain(int* punstable){/*{{{*/

	IssmDouble  penetration;
	bool        activate;
	int         unstable;
	IssmDouble  vx1;
	IssmDouble  vy1;
	IssmDouble  vx2;
	IssmDouble  vy2;

	/*Objects: */
	Tria  *tria1 = NULL;
	Tria  *tria2 = NULL;

	/*enum of element? */
	if(elements[0]->ObjectEnum()!=TriaEnum)_error_("only Tria element allowed for Riftfront load!");

	/*recover elements on both side of rift: */
	tria1=(Tria*)elements[0];
	tria2=(Tria*)elements[1];

	/*Is this constraint frozen? In which case we don't touch: */
	if (this->frozen){
		*punstable=0;
		return 1;
	}

	/*Is this rift segment state specified by user input? :*/
	if (this->state==OpenEnum || this->state==ClosedEnum){

		if(this->state==OpenEnum)this->active=0;
		if(this->state==ClosedEnum)this->active=1;

		/*this segment is like frozen, no instability here: */
		*punstable=0;
		return 1;
	}

	/*First recover velocity: */
	tria1->GetInputValue(&vx1,nodes[0],VxEnum);
	tria2->GetInputValue(&vx2,nodes[1],VxEnum);
	tria1->GetInputValue(&vy1,nodes[0],VyEnum);
	tria2->GetInputValue(&vy2,nodes[1],VyEnum);

	/*Node 1 faces node 2, compute penetration of 2 into 1 (V2-V1).N (with N normal vector, and V velocity vector: */
	penetration=(vx2-vx1)*normal[0]+(vy2-vy1)*normal[1];

	/*activation: */
	if(penetration<0)activate=true;
	else  activate=false;

	/*Here, we try to avoid zigzaging. When a penalty activates and deactivates for more than penalty_lock times, 
	 * we increase the fraction of melange:*/
	if(this->counter>this->penalty_lock){
		/*reset counter: */
		this->counter=0;
		/*increase melange fraction: */
		this->fraction+=fractionincrement;
		if (this->fraction>1)this->fraction=1.;
		//_printf_("riftfront " << this->Id() << " fraction: " << this->fraction << "\n");
	}

	//Figure out stability of this penalty
	if(this->active==activate){
		unstable=0;
	}
	else{
		unstable=1;
		this->counter++;
	}

	//Set penalty flag
	this->active=activate;

	//if ((penetration>0) && (this->active==1))_printf_("Riftfront " << Id() << " wants to be released\n");

	/*assign output pointer: */
	*punstable=unstable;
	return 1;
}
/*}}}*/
void   Riftfront::FreezeConstraints(void){/*{{{*/

	/*Just set frozen flag to 1: */
	this->frozen=1;

}
/*}}}*/
bool   Riftfront::IsFrozen(void){/*{{{*/

	/*Just set frozen flag to 1: */
	if(this->frozen)return 1;
	else return 0;
}
/*}}}*/
