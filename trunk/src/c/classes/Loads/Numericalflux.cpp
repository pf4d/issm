/*!\file Numericalflux.c
 * \brief: implementation of the Numericalflux object
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
#define NUMNODES_INTERNAL 4
#define NUMNODES_BOUNDARY 2

/*Numericalflux constructors and destructor*/
Numericalflux::Numericalflux(){/*{{{*/
	this->parameters = NULL;
	this->helement   = NULL;
	this->element    = NULL;
	this->hnodes     = NULL;
	this->hvertices  = NULL;
	this->nodes      = NULL;
}
/*}}}*/
Numericalflux::Numericalflux(int numericalflux_id,int i,int index,IoModel* iomodel, int in_analysis_type){/*{{{*/

	/* Intermediary */
	int  j;
	int  pos1,pos2,pos3,pos4;
	int  num_nodes;

	/*numericalflux constructor data: */
	int   numericalflux_elem_ids[2];
	int   numericalflux_vertex_ids[2];
	int   numericalflux_node_ids[4];
	int   numericalflux_type;

	/*Get edge*/
	int i1 = iomodel->faces[4*index+0];
	int i2 = iomodel->faces[4*index+1];
	int e1 = iomodel->faces[4*index+2];
	int e2 = iomodel->faces[4*index+3];

	/*First, see wether this is an internal or boundary edge (if e2=-1)*/
	if(e2==-1){
		/* Boundary edge, only one element */
		num_nodes=2;
		numericalflux_type=BoundaryEnum;
		numericalflux_elem_ids[0]=e1;
	}
	else{
		/* internal edge: connected to 2 elements */
		 num_nodes=4;
		numericalflux_type=InternalEnum;
		numericalflux_elem_ids[0]=e1;
		numericalflux_elem_ids[1]=e2;
	}

	/*1: Get vertices ids*/
	numericalflux_vertex_ids[0]=i1;
	numericalflux_vertex_ids[1]=i2;

	/*2: Get node ids*/
	if (numericalflux_type==InternalEnum){

		/*Now, we must get the nodes of the 4 nodes located on the edge*/

		/*2: Get the column where these ids are located in the index*/
		pos1=pos2=pos3=pos4=UNDEF;
		for(j=0;j<3;j++){
			if(iomodel->elements[3*(e1-1)+j]==i1) pos1=j+1;
			if(iomodel->elements[3*(e1-1)+j]==i2) pos2=j+1;
			if(iomodel->elements[3*(e2-1)+j]==i1) pos3=j+1;
			if(iomodel->elements[3*(e2-1)+j]==i2) pos4=j+1;
		}
		_assert_(pos1!=UNDEF && pos2!=UNDEF && pos3!=UNDEF && pos4!=UNDEF);

		/*3: We have the id of the elements and the position of the vertices in the index
		 * we can compute their dofs!*/
		numericalflux_node_ids[0]=iomodel->nodecounter+3*(e1-1)+pos1;
		numericalflux_node_ids[1]=iomodel->nodecounter+3*(e1-1)+pos2;
		numericalflux_node_ids[2]=iomodel->nodecounter+3*(e2-1)+pos3;
		numericalflux_node_ids[3]=iomodel->nodecounter+3*(e2-1)+pos4;
	}
	else{

		/*2: Get the column where these ids are located in the index*/
		pos1=pos2=UNDEF;
		for(j=0;j<3;j++){
			if(iomodel->elements[3*(e1-1)+j]==i1) pos1=j+1;
			if(iomodel->elements[3*(e1-1)+j]==i2) pos2=j+1;
		}
		_assert_(pos1!=UNDEF && pos2!=UNDEF);

		/*3: We have the id of the elements and the position of the vertices in the index
		 * we can compute their dofs!*/
		numericalflux_node_ids[0]=iomodel->nodecounter+3*(e1-1)+pos1;
		numericalflux_node_ids[1]=iomodel->nodecounter+3*(e1-1)+pos2;
	}

	/*Ok, we have everything to build the object: */
	this->id=numericalflux_id;
	this->analysis_type=in_analysis_type;
	this->flux_type = numericalflux_type;

	/*Hooks: */
	this->hnodes    =new Hook(numericalflux_node_ids,num_nodes);
	this->hvertices =new Hook(&numericalflux_vertex_ids[0],2);
	this->helement  =new Hook(numericalflux_elem_ids,1); // take only the first element for now

	//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
	this->parameters=NULL;
	this->element=NULL;
	this->nodes=NULL;
}
/*}}}*/
Numericalflux::~Numericalflux(){/*{{{*/
	this->parameters=NULL;
	delete helement;
	delete hnodes;
	delete hvertices;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Numericalflux::copy() {/*{{{*/

	Numericalflux* numericalflux=NULL;

	numericalflux=new Numericalflux();

	/*copy fields: */
	numericalflux->id=this->id;
	numericalflux->analysis_type=this->analysis_type;
	numericalflux->flux_type=this->flux_type;

	/*point parameters: */
	numericalflux->parameters=this->parameters;

	/*now deal with hooks and objects: */
	numericalflux->hnodes    = (Hook*)this->hnodes->copy();
	numericalflux->hvertices = (Hook*)this->hvertices->copy();
	numericalflux->helement  = (Hook*)this->helement->copy();

	/*corresponding fields*/
	numericalflux->nodes    = (Node**)numericalflux->hnodes->deliverp();
	numericalflux->vertices = (Vertex**)numericalflux->hvertices->deliverp();
	numericalflux->element  = (Element*)numericalflux->helement->delivers();

	return numericalflux;
}
/*}}}*/
void    Numericalflux::DeepEcho(void){/*{{{*/

	_printf_("Numericalflux:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	_printf_("   flux_type: " << this->flux_type<< "\n");
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
void    Numericalflux::Echo(void){/*{{{*/
	_printf_("Numericalflux:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	_printf_("   flux_type: " << this->flux_type<< "\n");
	hnodes->Echo();
	hvertices->Echo();
	helement->Echo();
	_printf_("   parameters: " << parameters << "\n");
}
/*}}}*/
int     Numericalflux::Id(void){/*{{{*/
	return id;
}
/*}}}*/
void    Numericalflux::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	_assert_(this);

	/*ok, marshall operations: */
	MARSHALLING_ENUM(NumericalfluxEnum);
	MARSHALLING(id);
	MARSHALLING(analysis_type);
	MARSHALLING(flux_type);

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
int     Numericalflux::ObjectEnum(void){/*{{{*/

	return NumericalfluxEnum;

}
/*}}}*/

/*Load virtual functions definitions:*/
void  Numericalflux::Configure(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

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
void  Numericalflux::CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/

	/*recover some parameters*/
	ElementMatrix* Ke=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	/*Just branch to the correct element stiffness matrix generator, according to the type of analysis we are carrying out: */
	switch(analysis_type){
		case MasstransportAnalysisEnum:
			Ke=CreateKMatrixMasstransport();
			break;
		case BalancethicknessAnalysisEnum:
			Ke=CreateKMatrixBalancethickness();
			break;
		case AdjointBalancethicknessAnalysisEnum:
			Ke=CreateKMatrixAdjointBalancethickness();
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
void  Numericalflux::CreatePVector(Vector<IssmDouble>* pf){/*{{{*/

	/*recover some parameters*/
	ElementVector* pe=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case MasstransportAnalysisEnum:
			pe=CreatePVectorMasstransport();
			break;
		case BalancethicknessAnalysisEnum:
			pe=CreatePVectorBalancethickness();
			break;
		case AdjointBalancethicknessAnalysisEnum:
			pe=CreatePVectorAdjointBalancethickness();
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
void  Numericalflux::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(nodes);

	switch(this->flux_type){
		case InternalEnum:
			for(int i=0;i<NUMNODES_INTERNAL;i++) lidlist[i]=nodes[i]->Lid();
			return;
		case BoundaryEnum:
			for(int i=0;i<NUMNODES_BOUNDARY;i++) lidlist[i]=nodes[i]->Lid();
			return;
		default:
			_error_("Numericalflux type " << EnumToStringx(this->flux_type) << " not supported yet");
	}
}
/*}}}*/
void  Numericalflux::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(nodes);

	switch(this->flux_type){
		case InternalEnum:
			for(int i=0;i<NUMNODES_INTERNAL;i++) sidlist[i]=nodes[i]->Sid();
			return;
		case BoundaryEnum:
			for(int i=0;i<NUMNODES_BOUNDARY;i++) sidlist[i]=nodes[i]->Sid();
			return;
		default:
			_error_("Numericalflux type " << EnumToStringx(this->flux_type) << " not supported yet");
	}
}
/*}}}*/
int   Numericalflux::GetNumberOfNodes(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return NUMNODES_INTERNAL;
		case BoundaryEnum:
			return NUMNODES_BOUNDARY;
		default:
			_error_("Numericalflux type " << EnumToStringx(this->flux_type) << " not supported yet");
	}

}
/*}}}*/
bool  Numericalflux::InAnalysis(int in_analysis_type){/*{{{*/
	if (in_analysis_type==this->analysis_type) return true;
	else return false;
}
/*}}}*/
bool  Numericalflux::IsPenalty(void){/*{{{*/
	return false;
}
/*}}}*/
void  Numericalflux::PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs,IssmDouble kmax){/*{{{*/

	/*No stiffness loads applied, do nothing: */
	return;

}
/*}}}*/
void  Numericalflux::PenaltyCreatePVector(Vector<IssmDouble>* pf,IssmDouble kmax){/*{{{*/

	/*No penalty loads applied, do nothing: */
	return;

}
/*}}}*/
void  Numericalflux::ResetHooks(){/*{{{*/

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
void  Numericalflux::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Numericalflux::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int set1_enum,int set2_enum){/*{{{*/

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

/*Numericalflux management*/
ElementMatrix* Numericalflux::CreateKMatrixAdjointBalancethickness(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return CreateKMatrixAdjointBalancethicknessInternal();
		case BoundaryEnum:
			return CreateKMatrixAdjointBalancethicknessBoundary();
		default:
			_error_("type not supported yet");
	}
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixAdjointBalancethicknessBoundary(void){/*{{{*/

	ElementMatrix* Ke=CreateKMatrixBalancethicknessBoundary();
	if(Ke) Ke->Transpose();
	return Ke;
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixAdjointBalancethicknessInternal(void){/*{{{*/

	ElementMatrix* Ke=CreateKMatrixBalancethicknessInternal();
	if (Ke) Ke->Transpose();
	return Ke;
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixBalancethickness(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return CreateKMatrixBalancethicknessInternal();
		case BoundaryEnum:
			return CreateKMatrixBalancethicknessBoundary();
		default:
			_error_("type not supported yet");
	}
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixBalancethicknessBoundary(void){/*{{{*/

	/* constants*/
	const int numdof=NDOF1*NUMNODES_BOUNDARY;

	/* Intermediaries*/
	int        i,j,ig,index1,index2;
	IssmDouble     DL,Jdet,vx,vy,mean_vx,mean_vy,UdotN;
	IssmDouble     xyz_list[NUMVERTICES][3];
	IssmDouble     normal[2];
	IssmDouble     L[numdof];
	IssmDouble     Ke_g[numdof][numdof];
	GaussTria *gauss;

	/*Initialize Element matrix and return if necessary*/
	ElementMatrix* Ke = NULL;
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* vxaverage_input=tria->inputs->GetInput(VxEnum);
	Input* vyaverage_input=tria->inputs->GetInput(VyEnum);
	GetNormal(&normal[0],xyz_list);

	/*Check wether it is an inflow or outflow BC (0 is the middle of the segment)*/
	index1=tria->GetNodeIndex(nodes[0]);
	index2=tria->GetNodeIndex(nodes[1]);

	gauss=new GaussTria();
	gauss->GaussEdgeCenter(index1,index2);
	vxaverage_input->GetInputValue(&mean_vx,gauss);
	vyaverage_input->GetInputValue(&mean_vy,gauss);
	delete gauss;

	UdotN=mean_vx*normal[0]+mean_vy*normal[1];
	if (UdotN<=0){
		return NULL; /*(u,n)<0 -> inflow, PenaltyCreatePVector will take care of it*/
	}
	else{
		Ke=new ElementMatrix(nodes,NUMNODES_BOUNDARY,this->parameters);
	}

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(index1,index2,2);
	for(ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		tria->GetSegmentNodalFunctions(&L[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		DL=gauss->weight*Jdet*UdotN;

		TripleMultiply(&L[0],1,numdof,1,
					&DL,1,1,0,
					&L[0],1,numdof,0,
					&Ke_g[0][0],0);

		for(i=0;i<numdof;i++) for(j=0;j<numdof;j++) Ke->values[i*numdof+j]+=Ke_g[i][j];
	} 

	/*Clean up and return*/
	delete gauss;
	return Ke;
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixBalancethicknessInternal(void){/*{{{*/

	/* constants*/
	const int numdof=NDOF1*NUMNODES_INTERNAL;

	/* Intermediaries*/
	int        i,j,ig,index1,index2;
	IssmDouble     DL1,DL2,Jdet,vx,vy,UdotN;
	IssmDouble     xyz_list[NUMVERTICES][3];
	IssmDouble     normal[2];
	IssmDouble     B[numdof];
	IssmDouble     Bprime[numdof];
	IssmDouble     Ke_g1[numdof][numdof];
	IssmDouble     Ke_g2[numdof][numdof];
	GaussTria *gauss;

	/*Initialize Element matrix and return if necessary*/
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;
	ElementMatrix* Ke=new ElementMatrix(nodes,NUMNODES_INTERNAL,this->parameters);

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* vxaverage_input=tria->inputs->GetInput(VxEnum);
	Input* vyaverage_input=tria->inputs->GetInput(VyEnum);
	GetNormal(&normal[0],xyz_list);

	/* Start  looping on the number of gaussian points: */
	index1=tria->GetNodeIndex(nodes[0]);
	index2=tria->GetNodeIndex(nodes[1]);
	gauss=new GaussTria(index1,index2,2);
	for(ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		tria->GetSegmentBFlux(&B[0],gauss,index1,index2,tria->FiniteElement());
		tria->GetSegmentBprimeFlux(&Bprime[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		DL1=gauss->weight*Jdet*UdotN/2;
		DL2=gauss->weight*Jdet*fabs(UdotN)/2;

		TripleMultiply(&B[0],1,numdof,1,
					&DL1,1,1,0,
					&Bprime[0],1,numdof,0,
					&Ke_g1[0][0],0);
		TripleMultiply(&B[0],1,numdof,1,
					&DL2,1,1,0,
					&B[0],1,numdof,0,
					&Ke_g2[0][0],0);

		for(i=0;i<numdof;i++) for(j=0;j<numdof;j++) Ke->values[i*numdof+j]+=Ke_g1[i][j];
		for(i=0;i<numdof;i++) for(j=0;j<numdof;j++) Ke->values[i*numdof+j]+=Ke_g2[i][j];
	}

	/*Clean up and return*/
	delete gauss;
	return Ke;
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixMasstransport(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return CreateKMatrixMasstransportInternal();
		case BoundaryEnum:
			return CreateKMatrixMasstransportBoundary();
		default:
			_error_("type not supported yet");
	}
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixMasstransportBoundary(void){/*{{{*/

	/* constants*/
	const int numdof=NDOF1*NUMNODES_BOUNDARY;

	/* Intermediaries*/
	int        i,j,ig,index1,index2;
	IssmDouble     DL,Jdet,dt,vx,vy,mean_vx,mean_vy,UdotN;
	IssmDouble     xyz_list[NUMVERTICES][3];
	IssmDouble     normal[2];
	IssmDouble     L[numdof];
	IssmDouble     Ke_g[numdof][numdof];
	GaussTria *gauss;

	/*Initialize Element matrix and return if necessary*/
	ElementMatrix* Ke = NULL;
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	Input* vxaverage_input=tria->inputs->GetInput(VxEnum); _assert_(vxaverage_input);
	Input* vyaverage_input=tria->inputs->GetInput(VyEnum); _assert_(vyaverage_input);
	GetNormal(&normal[0],xyz_list);

	/*Check wether it is an inflow or outflow BC (0 is the middle of the segment)*/
	index1=tria->GetNodeIndex(nodes[0]);
	index2=tria->GetNodeIndex(nodes[1]);

	gauss=new GaussTria();
	gauss->GaussEdgeCenter(index1,index2);
	vxaverage_input->GetInputValue(&mean_vx,gauss);
	vyaverage_input->GetInputValue(&mean_vy,gauss);
	delete gauss;

	UdotN=mean_vx*normal[0]+mean_vy*normal[1];
	if (UdotN<=0){
		return NULL; /*(u,n)<0 -> inflow, PenaltyCreatePVector will take care of it*/
	}
	else{
		Ke=new ElementMatrix(nodes,NUMNODES_BOUNDARY,this->parameters);
	}

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(index1,index2,2);
	for(ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		tria->GetSegmentNodalFunctions(&L[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		DL=gauss->weight*Jdet*dt*UdotN;

		TripleMultiply(&L[0],1,numdof,1,
					&DL,1,1,0,
					&L[0],1,numdof,0,
					&Ke_g[0][0],0);

		for(i=0;i<numdof;i++) for(j=0;j<numdof;j++) Ke->values[i*numdof+j]+=Ke_g[i][j];
	} 

	/*Clean up and return*/
	delete gauss;
	return Ke;
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixMasstransportInternal(void){/*{{{*/

	/* constants*/
	const int numdof=NDOF1*NUMNODES_INTERNAL;

	/* Intermediaries*/
	int        i,j,ig,index1,index2;
	IssmDouble     DL1,DL2,Jdet,dt,vx,vy,UdotN;
	IssmDouble     xyz_list[NUMVERTICES][3];
	IssmDouble     normal[2];
	IssmDouble     B[numdof];
	IssmDouble     Bprime[numdof];
	IssmDouble     Ke_g1[numdof][numdof];
	IssmDouble     Ke_g2[numdof][numdof];
	GaussTria *gauss;

	/*Initialize Element matrix and return if necessary*/
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;
	ElementMatrix* Ke=new ElementMatrix(nodes,NUMNODES_INTERNAL,this->parameters);

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	Input* vxaverage_input=tria->inputs->GetInput(VxEnum);
	Input* vyaverage_input=tria->inputs->GetInput(VyEnum);
	GetNormal(&normal[0],xyz_list);

	/* Start  looping on the number of gaussian points: */
	index1=tria->GetNodeIndex(nodes[0]);
	index2=tria->GetNodeIndex(nodes[1]);
	gauss=new GaussTria(index1,index2,2);
	for(ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		tria->GetSegmentBFlux(&B[0],gauss,index1,index2,tria->FiniteElement());
		tria->GetSegmentBprimeFlux(&Bprime[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		DL1=gauss->weight*Jdet*dt*UdotN/2;
		DL2=gauss->weight*Jdet*dt*fabs(UdotN)/2;

		TripleMultiply(&B[0],1,numdof,1,
					&DL1,1,1,0,
					&Bprime[0],1,numdof,0,
					&Ke_g1[0][0],0);
		TripleMultiply(&B[0],1,numdof,1,
					&DL2,1,1,0,
					&B[0],1,numdof,0,
					&Ke_g2[0][0],0);

		for(i=0;i<numdof;i++) for(j=0;j<numdof;j++) Ke->values[i*numdof+j]+=Ke_g1[i][j];
		for(i=0;i<numdof;i++) for(j=0;j<numdof;j++) Ke->values[i*numdof+j]+=Ke_g2[i][j];
	}

	/*Clean up and return*/
	delete gauss;
	return Ke;
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorAdjointBalancethickness(void){/*{{{*/

	/*No PVector for the Adjoint*/
	return NULL;
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorBalancethickness(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return CreatePVectorBalancethicknessInternal();
		case BoundaryEnum:
			return CreatePVectorBalancethicknessBoundary();
		default:
			_error_("type not supported yet");
	}
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorBalancethicknessBoundary(void){/*{{{*/

	/* constants*/
	const int numdof=NDOF1*NUMNODES_BOUNDARY;

	/* Intermediaries*/
	int        i,ig,index1,index2;
	IssmDouble DL,Jdet,vx,vy,mean_vx,mean_vy,UdotN,thickness;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble normal[2];
	IssmDouble L[numdof];
	GaussTria *gauss;

	/*Initialize Load Vector and return if necessary*/
	ElementVector* pe = NULL;
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* vxaverage_input=tria->inputs->GetInput(VxEnum); _assert_(vxaverage_input); 
	Input* vyaverage_input=tria->inputs->GetInput(VyEnum); _assert_(vyaverage_input);
	Input* thickness_input=tria->inputs->GetInput(ThicknessEnum); _assert_(thickness_input);
	GetNormal(&normal[0],xyz_list);

	/*Check wether it is an inflow or outflow BC (0 is the middle of the segment)*/
	index1=tria->GetNodeIndex(nodes[0]);
	index2=tria->GetNodeIndex(nodes[1]);

	gauss=new GaussTria();
	gauss->GaussEdgeCenter(index1,index2);
	vxaverage_input->GetInputValue(&mean_vx,gauss);
	vyaverage_input->GetInputValue(&mean_vy,gauss);
	delete gauss;
	UdotN=mean_vx*normal[0]+mean_vy*normal[1];
	if (UdotN>0){
		return NULL; /*(u,n)>0 -> outflow, PenaltyCreateKMatrix will take care of it*/
	}
	else{
		pe=new ElementVector(nodes,NUMNODES_BOUNDARY,this->parameters);
	}

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(index1,index2,2);
	for(ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		tria->GetSegmentNodalFunctions(&L[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		thickness_input->GetInputValue(&thickness,gauss);

		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		DL= - gauss->weight*Jdet*UdotN*thickness;

		for(i=0;i<numdof;i++) pe->values[i] += DL*L[i];
	}

	/*Clean up and return*/
	delete gauss;
	return pe;
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorBalancethicknessInternal(void){/*{{{*/

	/*Nothing added to PVector*/
	return NULL;

}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorMasstransport(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return CreatePVectorMasstransportInternal();
		case BoundaryEnum:
			return CreatePVectorMasstransportBoundary();
		default:
			_error_("type not supported yet");
	}
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorMasstransportBoundary(void){/*{{{*/

	/* constants*/
	const int numdof=NDOF1*NUMNODES_BOUNDARY;

	/* Intermediaries*/
	int        i,ig,index1,index2;
	IssmDouble     DL,Jdet,dt,vx,vy,mean_vx,mean_vy,UdotN,thickness;
	IssmDouble     xyz_list[NUMVERTICES][3];
	IssmDouble     normal[2];
	IssmDouble     L[numdof];
	GaussTria *gauss;

	/*Initialize Load Vector and return if necessary*/
	ElementVector* pe = NULL;
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	Input* vxaverage_input   =tria->inputs->GetInput(VxEnum);                     _assert_(vxaverage_input); 
	Input* vyaverage_input   =tria->inputs->GetInput(VyEnum);                     _assert_(vyaverage_input);
	Input* spcthickness_input=tria->inputs->GetInput(MasstransportSpcthicknessEnum); _assert_(spcthickness_input);
	GetNormal(&normal[0],xyz_list);

	/*Check wether it is an inflow or outflow BC (0 is the middle of the segment)*/
	index1=tria->GetNodeIndex(nodes[0]);
	index2=tria->GetNodeIndex(nodes[1]);

	gauss=new GaussTria();
	gauss->GaussEdgeCenter(index1,index2);
	vxaverage_input->GetInputValue(&mean_vx,gauss);
	vyaverage_input->GetInputValue(&mean_vy,gauss);
	delete gauss;

	UdotN=mean_vx*normal[0]+mean_vy*normal[1];
	if (UdotN>0){
		return NULL; /*(u,n)>0 -> outflow, PenaltyCreateKMatrix will take care of it*/
	}
	else{
		pe=new ElementVector(nodes,NUMNODES_BOUNDARY,this->parameters);
	}

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(index1,index2,2);
	for(ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		tria->GetSegmentNodalFunctions(&L[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		spcthickness_input->GetInputValue(&thickness,gauss);
		if(xIsNan<IssmDouble>(thickness)) _error_("Cannot weakly apply constraint because NaN was provided");

		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		DL= - gauss->weight*Jdet*dt*UdotN*thickness;

		for(i=0;i<numdof;i++) pe->values[i] += DL*L[i];
	}

	/*Clean up and return*/
	delete gauss;
	return pe;
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorMasstransportInternal(void){/*{{{*/

	/*Nothing added to PVector*/
	return NULL;

}
/*}}}*/
void           Numericalflux::GetNormal(IssmDouble* normal,IssmDouble xyz_list[4][3]){/*{{{*/

	/*Build unit outward pointing vector*/
	IssmDouble vector[2];
	IssmDouble norm;

	vector[0]=xyz_list[1][0] - xyz_list[0][0];
	vector[1]=xyz_list[1][1] - xyz_list[0][1];

	norm=sqrt(pow(vector[0],2.0)+pow(vector[1],2.0));

	normal[0]= + vector[1]/norm;
	normal[1]= - vector[0]/norm;
}
/*}}}*/
