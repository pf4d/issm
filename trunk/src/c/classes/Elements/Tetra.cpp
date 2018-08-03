/*!\file Tetra.cpp
 * \brief: implementation of the Tetrament object
 */
/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>
#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*Element macros*/
#define NUMVERTICES 4

/*Constructors/destructor/copy*/
Tetra::Tetra(int seg_id, int seg_sid, int index, IoModel* iomodel,int nummodels)/*{{{*/
		:ElementHook(nummodels,index+1,NUMVERTICES,iomodel){

			/*id: */
			this->id  = seg_id;
			this->sid = seg_sid;

			//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
			this->parameters = NULL;

			/*intialize inputs: */
			this->inputs  = new Inputs();

			/*initialize pointers:*/
			this->nodes    = NULL;
			this->vertices = NULL;
			this->material = NULL;
			this->matpar   = NULL;

			/*Only allocate pointer*/
			this->element_type_list=xNew<int>(nummodels);
		}
/*}}}*/
Tetra::~Tetra(){/*{{{*/
	this->parameters=NULL;
}
/*}}}*/
Object* Tetra::copy() {/*{{{*/

	int i;
	Tetra* tetra=NULL;

	tetra=new Tetra();

	//deal with TetraRef mother class
	int nanalyses = this->numanalyses;
	if(nanalyses > 0){
		tetra->element_type_list=xNew<int>(nanalyses);
		for(i=0;i<nanalyses;i++){
			if (this->element_type_list[i]) tetra->element_type_list[i]=this->element_type_list[i];
			else tetra->element_type_list[i] = 0;
		}
	}
	else tetra->element_type_list = NULL;
	tetra->element_type=this->element_type;
	tetra->numanalyses=nanalyses;

	//deal with ElementHook mother class
	if (this->hnodes){
		tetra->hnodes=xNew<Hook*>(tetra->numanalyses);
		for(i=0;i<tetra->numanalyses;i++){
			if (this->hnodes[i]) tetra->hnodes[i] = (Hook*)(this->hnodes[i]->copy());
			else tetra->hnodes[i] = NULL;
		}
	}
	else tetra->hnodes = NULL;

	tetra->hvertices = (Hook*)this->hvertices->copy();
	tetra->hmaterial = (Hook*)this->hmaterial->copy();
	tetra->hmatpar   = (Hook*)this->hmatpar->copy();
	tetra->hneighbors = NULL;

	/*deal with Tria fields: */
	tetra->id  = this->id;
	tetra->sid = this->sid;
	if(this->inputs) tetra->inputs = (Inputs*)(this->inputs->Copy());
	else tetra->inputs=new Inputs();

	/*point parameters: */
	tetra->parameters=this->parameters;

	/*recover objects: */
	unsigned int num_nodes = 4;
	tetra->nodes = xNew<Node*>(num_nodes); //we cannot rely on an analysis_counter to tell us which analysis_type we are running, so we just copy the nodes.
	for(i=0;i<num_nodes;i++) if(this->nodes[i]) tetra->nodes[i]=this->nodes[i]; else tetra->nodes[i] = NULL;

	tetra->vertices = (Vertex**)this->hvertices->deliverp();
	tetra->material = (Material*)this->hmaterial->delivers();
	tetra->matpar   = (Matpar*)this->hmatpar->delivers();

	return tetra;
}
/*}}}*/
void Tetra::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(TetraEnum);

	/*Call parent classes: */
	ElementHook::Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	Element::MarshallElement(pmarshalled_data,pmarshalled_data_size,marshall_direction,this->numanalyses);
	TetraRef::Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);

	vertices = (Vertex**)this->hvertices->deliverp();
	material = (Material*)this->hmaterial->delivers();
	matpar   = (Matpar*)this->hmatpar->delivers();

}
/*}}}*/

void     Tetra::AddInput(int input_enum,IssmDouble* values, int interpolation_enum){/*{{{*/

	/*Call inputs method*/
	_assert_(this->inputs);
	this->inputs->AddInput(new TetraInput(input_enum,values,interpolation_enum));
}
/*}}}*/
void     Tetra::Configure(Elements* elementsin, Loads* loadsin, Nodes* nodesin,Vertices* verticesin, Materials* materialsin, Parameters* parametersin){/*{{{*/

	int analysis_counter;

	/*go into parameters and get the analysis_counter: */
	parametersin->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Get Element type*/
	this->element_type=this->element_type_list[analysis_counter];

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective 
	 * datasets, using internal ids and offsets hidden in hooks: */
	if (this->hnodes[analysis_counter]) this->hnodes[analysis_counter]->configure(nodesin);
	this->hvertices->configure(verticesin);
	this->hmaterial->configure(materialsin);
	this->hmatpar->configure(materialsin);

	/*Now, go pick up the objects inside the hooks: */
	if (this->hnodes[analysis_counter]) this->nodes=(Node**)this->hnodes[analysis_counter]->deliverp();
	else this->nodes=NULL;
	this->vertices          = (Vertex**)this->hvertices->deliverp();
	this->material          = (Material*)this->hmaterial->delivers();
	this->matpar            = (Matpar*)this->hmatpar->delivers();

	/*point parameters to real dataset: */
	this->parameters=parametersin;

	/*get inputs configured too: */
	this->inputs->Configure(parameters);
}
/*}}}*/
void     Tetra::ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz){/*{{{*/

	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble xmin,ymin,zmin;
	IssmDouble xmax,ymax,zmax;

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	xmin=xyz_list[0][0]; xmax=xyz_list[0][0];
	ymin=xyz_list[0][1]; ymax=xyz_list[0][1];
	zmin=xyz_list[0][2]; zmax=xyz_list[0][2];

	for(int i=1;i<NUMVERTICES;i++){
		if(xyz_list[i][0]<xmin) xmin=xyz_list[i][0];
		if(xyz_list[i][0]>xmax) xmax=xyz_list[i][0];
		if(xyz_list[i][1]<ymin) ymin=xyz_list[i][1];
		if(xyz_list[i][1]>ymax) ymax=xyz_list[i][1];
		if(xyz_list[i][2]<zmin) zmin=xyz_list[i][2];
		if(xyz_list[i][2]>zmax) zmax=xyz_list[i][2];
	}

	*hx=xmax-xmin;
	*hy=ymax-ymin;
	*hz=zmax-zmin;
}
/*}}}*/
void     Tetra::FaceOnBaseIndices(int* pindex1,int* pindex2,int* pindex3){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[4][3] = {{0,1,2},{0,3,1},{1,3,2},{0,2,3}};

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);

	for(int i=0;i<4;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1. && values[indices[i][2]] == 1.){
			*pindex1 = indices[i][0];
			*pindex2 = indices[i][1];
			*pindex3 = indices[i][2];
			return;
		}
	}

	_error_("Could not find 3 vertices on bed");
}
/*}}}*/
void     Tetra::FaceOnFrontIndices(int* pindex1,int* pindex2,int* pindex3){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[4][3] = {{0,1,2},{0,3,1},{1,3,2},{0,2,3}};

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MaskIceLevelsetEnum);

	for(int i=0;i<4;i++){
		if(values[indices[i][0]] == 0. && values[indices[i][1]] == 0. && values[indices[i][2]] == 0.){
			*pindex1 = indices[i][0];
			*pindex2 = indices[i][1];
			*pindex3 = indices[i][2];
			return;
		}
	}

	_error_("Could not find 3 vertices on bed");
}
/*}}}*/
void     Tetra::FaceOnSurfaceIndices(int* pindex1,int* pindex2,int* pindex3){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[4][3] = {{0,1,2},{0,3,1},{1,3,2},{0,2,3}};

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);

	for(int i=0;i<4;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1. && values[indices[i][2]] == 1.){
			*pindex1 = indices[i][0];
			*pindex2 = indices[i][1];
			*pindex3 = indices[i][2];
			return;
		}
	}

	_error_("Could not find 3 vertices on bed");
}
/*}}}*/
int      Tetra::FiniteElement(void){/*{{{*/
	return this->element_type;
} /*}}}*/
int      Tetra::GetElementType(){/*{{{*/

	/*return TetraRef field*/
	return this->element_type;
}
/*}}}*/
void     Tetra::GetInputValue(IssmDouble* pvalue,Node* node,int enumtype){/*{{{*/

	Input* input=inputs->GetInput(enumtype);
	if(!input) _error_("No input of type " << EnumToStringx(enumtype) << " found in tria");

	GaussTetra* gauss=new GaussTetra();
	gauss->GaussVertex(this->GetNodeIndex(node));

	input->GetInputValue(pvalue,gauss);
	delete gauss;
}
/*}}}*/
int      Tetra::GetNodeIndex(Node* node){/*{{{*/

	_assert_(nodes);
	int numnodes = this->NumberofNodes(this->element_type);

	for(int i=0;i<numnodes;i++){
		if(node==nodes[i]) return i;
	}
	_error_("Node provided not found among element nodes");

}
/*}}}*/
int      Tetra::GetNumberOfNodes(void){/*{{{*/
	return this->NumberofNodes(this->element_type);
}
/*}}}*/
int      Tetra::GetNumberOfVertices(void){/*{{{*/
	return NUMVERTICES;
}
/*}}}*/
void     Tetra::GetVerticesCoordinatesBase(IssmDouble** pxyz_list){/*{{{*/

	int        indices[3];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Element XYZ list*/
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);

	/*Allocate Output*/
	IssmDouble* xyz_list_edge = xNew<IssmDouble>(3*3);
	this->FaceOnBaseIndices(&indices[0],&indices[1],&indices[2]);
	for(int i=0;i<3;i++) for(int j=0;j<3;j++) xyz_list_edge[i*3+j]=xyz_list[indices[i]][j];

	/*Assign output pointer*/
	*pxyz_list = xyz_list_edge;

}/*}}}*/
void     Tetra::GetVerticesCoordinatesTop(IssmDouble** pxyz_list){/*{{{*/

	int        indices[3];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Element XYZ list*/
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);

	/*Allocate Output*/
	IssmDouble* xyz_list_edge = xNew<IssmDouble>(3*3);
	this->FaceOnSurfaceIndices(&indices[0],&indices[1],&indices[2]);
	for(int i=0;i<3;i++) for(int j=0;j<3;j++) xyz_list_edge[i*3+j]=xyz_list[indices[i]][j];

	/*Assign output pointer*/
	*pxyz_list = xyz_list_edge;

}/*}}}*/
bool     Tetra::HasFaceOnBase(){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);
	sum = values[0]+values[1]+values[2]+values[3];

	_assert_(sum==0. || sum==1. || sum==2. || sum==3.);

	if(sum==3){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
bool     Tetra::HasFaceOnSurface(){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);
	sum = values[0]+values[1]+values[2]+values[3];

	_assert_(sum==0. || sum==1. || sum==2. || sum==3.);

	if(sum==3){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
void     Tetra::InputUpdateFromIoModel(int index,IoModel* iomodel){ /*{{{*/

	/*Intermediaries*/
	int         i,j;
	int         tetra_vertex_ids[NUMVERTICES];
	IssmDouble  nodeinputs[NUMVERTICES];
	IssmDouble  cmmininputs[NUMVERTICES];
	IssmDouble  cmmaxinputs[NUMVERTICES];

	IssmDouble  yts;
	bool    control_analysis;
	char**  controls = NULL;
	int     num_control_type,num_responses;

	/*Fetch parameters: */
	iomodel->FindConstant(&yts,"md.constants.yts");
	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");
	if(control_analysis) iomodel->FindConstant(&num_control_type,"md.inversion.num_control_parameters");
	if(control_analysis) iomodel->FindConstant(&num_responses,"md.inversion.num_cost_functions");

	/*Recover vertices ids needed to initialize inputs*/
	_assert_(iomodel->elements);
	for(i=0;i<NUMVERTICES;i++){ 
		tetra_vertex_ids[i]=iomodel->elements[NUMVERTICES*index+i]; //ids for vertices are in the elements array from Matlab
	}

	/*Control Inputs*/
	if (control_analysis){
		iomodel->FindConstant(&controls,NULL,"md.inversion.control_parameters");
		for(i=0;i<num_control_type;i++){
			_assert_(controls[i]);
			int control = StringToEnumx(controls[i]);
			switch(control){
				case BalancethicknessThickeningRateEnum:
					if (iomodel->Data("md.balancethickness.thickening_rate")){
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.balancethickness.thickening_rate")[tetra_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i]/yts;
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i]/yts;
						this->inputs->AddInput(new ControlInput(BalancethicknessThickeningRateEnum,TetraInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case VxEnum:
					if (iomodel->Data("md.initialization.vx")){
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.initialization.vx")[tetra_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i]/yts;
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i]/yts;
						this->inputs->AddInput(new ControlInput(VxEnum,TetraInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case VyEnum:
					if (iomodel->Data("md.initialization.vy")){
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.initialization.vy")[tetra_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i]/yts;
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i]/yts;
						this->inputs->AddInput(new ControlInput(VyEnum,TetraInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case FrictionCoefficientEnum:
					if (iomodel->Data("md.friction.coefficient")){
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.friction.coefficient")[tetra_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i];
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i];
						this->inputs->AddInput(new ControlInput(FrictionCoefficientEnum,TetraInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case MaterialsRheologyBbarEnum:
					if(iomodel->Data("md.materials.rheology_B")){
						for(j=0;j<NUMVERTICES;j++) nodeinputs[j]=iomodel->Data("md.materials.rheology_B")[tetra_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i];
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i];
						this->inputs->AddInput(new ControlInput(MaterialsRheologyBEnum,TetraInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case DamageDbarEnum:
					if(iomodel->Data("md.damage.D")){
						for(j=0;j<NUMVERTICES;j++) nodeinputs[j]=iomodel->Data("md.damage.D")[tetra_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i];
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tetra_vertex_ids[j]-1)*num_control_type+i];
						this->inputs->AddInput(new ControlInput(DamageDEnum,TetraInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				default:
					_error_("Control " << EnumToStringx(control) << " not implemented yet");
			}
		}
		for(i=0;i<num_control_type;i++) xDelete<char>(controls[i]);
		xDelete<char*>(controls);
	}

	/*Need to know the type of approximation for this element*/
	if(iomodel->Data("md.flowequation.element_equation")){
		this->inputs->AddInput(new IntInput(ApproximationEnum,IoCodeToEnumElementEquation(reCast<int>(iomodel->Data("md.flowequation.element_equation")[index]))));
	}

	/*DatasetInputs*/
	if (control_analysis && iomodel->Data("md.inversion.cost_functions_coefficients")) {

		/*Generate cost functions associated with the iomodel*/
		char**	cost_functions			= NULL;
		int*		cost_functions_enums = NULL;
		int		num_cost_functions;

		iomodel->FindConstant(&num_cost_functions,"md.inversion.num_cost_functions");
		iomodel->FindConstant(&cost_functions,&num_cost_functions,"md.inversion.cost_functions");
		if(num_cost_functions<1) _error_("No cost functions found");
		cost_functions_enums=xNew<int>(num_cost_functions);
		for(j=0;j<num_cost_functions;j++){ cost_functions_enums[j]=StringToEnumx(cost_functions[j]); }

		/*Create inputs and add to DataSetInput*/
		DatasetInput* datasetinput=new DatasetInput(InversionCostFunctionsCoefficientsEnum);
		for(i=0;i<num_responses;i++){
			for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.inversion.cost_functions_coefficients")[(tetra_vertex_ids[j]-1)*num_responses+i];
			datasetinput->AddInput(new TetraInput(InversionCostFunctionsCoefficientsEnum,nodeinputs,P1Enum),cost_functions_enums[i]);
		}

		/*Add datasetinput to element inputs*/
		this->inputs->AddInput(datasetinput);

		/*Clean up cost functions*/
		xDelete<int>(cost_functions_enums);
		for(int j=0;j<num_cost_functions;j++) xDelete<char>(cost_functions[j]); 
		xDelete<char*>(cost_functions);
	}
}
/*}}}*/
void     Tetra::InputUpdateFromSolutionOneDof(IssmDouble* solution,int enum_type){/*{{{*/

	/*Intermediary*/
	int* doflist = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(this->element_type);

	/*Fetch dof list and allocate solution vector*/
	GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values    = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		values[i]=solution[doflist[i]];
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*Add input to the element: */
	this->inputs->AddInput(new TetraInput(enum_type,values,P1Enum));

	/*Free ressources:*/
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);
}
/*}}}*/
bool     Tetra::IsIcefront(void){/*{{{*/

	/*Retrieve all inputs and parameters*/
	IssmDouble ls[NUMVERTICES];
	GetInputListOnVertices(&ls[0],MaskIceLevelsetEnum);

	/* If only one vertex has ice, there is an ice front here */
	if(IsIceInElement()){
		int nrice=0;       
		for(int i=0;i<NUMVERTICES;i++) if(ls[i]<0.) nrice++;
		if(nrice==1) return true;
	}
	return false;
}/*}}}*/
bool     Tetra::IsOnBase(){/*{{{*/
	return HasFaceOnBase();
}
/*}}}*/
bool     Tetra::IsOnSurface(){/*{{{*/
	return HasFaceOnSurface();
}
/*}}}*/
void     Tetra::JacobianDeterminant(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetJacobianDeterminant(pJdet,xyz_list,(GaussTetra*)gauss);

}
/*}}}*/
void     Tetra::JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetJacobianDeterminantFace(pJdet,xyz_list_base,(GaussTetra*)gauss);

}
/*}}}*/
void     Tetra::JacobianDeterminantSurface(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetJacobianDeterminantFace(pJdet,xyz_list,(GaussTetra*)gauss);

}
/*}}}*/
void     Tetra::JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetJacobianDeterminantFace(pJdet,xyz_list_base,(GaussTetra*)gauss);

}
/*}}}*/
Gauss*   Tetra::NewGauss(void){/*{{{*/
	return new GaussTetra();
}
/*}}}*/
Gauss*   Tetra::NewGauss(int order){/*{{{*/
	return new GaussTetra(order);
}
/*}}}*/
Gauss*   Tetra::NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert){/*{{{*/
	/*FIXME: this is messed up, should provide indices and not xyz_list!*/
	int indices[3];
	this->FaceOnFrontIndices(&indices[0],&indices[1],&indices[2]);
	return new GaussTetra(indices[0],indices[1],indices[2],max(order_horiz,order_vert));
}
/*}}}*/
Gauss*   Tetra::NewGaussBase(int order){/*{{{*/

	int indices[3];
	this->FaceOnBaseIndices(&indices[0],&indices[1],&indices[2]);
	return new GaussTetra(indices[0],indices[1],indices[2],order);
}
/*}}}*/
Gauss*   Tetra::NewGaussTop(int order){/*{{{*/

	int indices[3];
	this->FaceOnSurfaceIndices(&indices[0],&indices[1],&indices[2]);
	return new GaussTetra(indices[0],indices[1],indices[2],order);
}
/*}}}*/
void     Tetra::NodalFunctions(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctions(basis,(GaussTetra*)gauss,this->element_type);

}
/*}}}*/
void     Tetra::NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussTetra*)gauss,this->element_type);

}
/*}}}*/
void     Tetra::NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussTetra*)gauss,this->VelocityInterpolation());

}
/*}}}*/
void     Tetra::NodalFunctionsPressure(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctions(basis,(GaussTetra*)gauss,this->PressureInterpolation());

}
/*}}}*/
void     Tetra::NodalFunctionsTensor(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctions(basis,(GaussTetra*)gauss,this->TensorInterpolation());

}
/*}}}*/
void     Tetra::NodalFunctionsVelocity(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctions(basis,(GaussTetra*)gauss,this->VelocityInterpolation());

}
/*}}}*/
void     Tetra::NormalBase(IssmDouble* bed_normal,IssmDouble* xyz_list){/*{{{*/

	IssmDouble v13[3],v23[3];
	IssmDouble normal[3];
	IssmDouble normal_norm;

	for(int i=0;i<3;i++){
		v13[i]=xyz_list[0*3+i]-xyz_list[2*3+i];
		v23[i]=xyz_list[1*3+i]-xyz_list[2*3+i];
	}

	normal[0]=v13[1]*v23[2]-v13[2]*v23[1];
	normal[1]=v13[2]*v23[0]-v13[0]*v23[2];
	normal[2]=v13[0]*v23[1]-v13[1]*v23[0];
	normal_norm=sqrt(normal[0]*normal[0]+ normal[1]*normal[1]+ normal[2]*normal[2]);

	/*Bed normal is opposite to surface normal*/
	bed_normal[0]=-normal[0]/normal_norm;
	bed_normal[1]=-normal[1]/normal_norm;
	bed_normal[2]=-normal[2]/normal_norm;

	_assert_(bed_normal[2]<0.);
}
/*}}}*/
void     Tetra::NormalSection(IssmDouble* normal,IssmDouble* xyz_list){/*{{{*/

	/*Build unit outward pointing vector*/
	IssmDouble AB[3];
	IssmDouble AC[3];
	IssmDouble norm;

	AB[0]=xyz_list[1*3+0] - xyz_list[0*3+0];
	AB[1]=xyz_list[1*3+1] - xyz_list[0*3+1];
	AB[2]=xyz_list[1*3+2] - xyz_list[0*3+2];
	AC[0]=xyz_list[2*3+0] - xyz_list[0*3+0];
	AC[1]=xyz_list[2*3+1] - xyz_list[0*3+1];
	AC[2]=xyz_list[2*3+2] - xyz_list[0*3+2];

	cross(normal,AB,AC);
	norm=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);

	for(int i=0;i<3;i++) normal[i]=normal[i]/norm;
}
/*}}}*/
void     Tetra::NormalTop(IssmDouble* top_normal,IssmDouble* xyz_list){/*{{{*/

	IssmDouble v13[3],v23[3];
	IssmDouble normal[3];
	IssmDouble normal_norm;

	for(int i=0;i<3;i++){
		v13[i]=xyz_list[0*3+i]-xyz_list[2*3+i];
		v23[i]=xyz_list[1*3+i]-xyz_list[2*3+i];
	}

	normal[0]=v13[1]*v23[2]-v13[2]*v23[1];
	normal[1]=v13[2]*v23[0]-v13[0]*v23[2];
	normal[2]=v13[0]*v23[1]-v13[1]*v23[0];
	normal_norm=sqrt(normal[0]*normal[0]+ normal[1]*normal[1]+ normal[2]*normal[2]);

	top_normal[0]=normal[0]/normal_norm;
	top_normal[1]=normal[1]/normal_norm;
	top_normal[2]=normal[2]/normal_norm;
	_assert_(top_normal[2]>0.);
}
/*}}}*/
int      Tetra::NumberofNodesPressure(void){/*{{{*/
	return TetraRef::NumberofNodes(this->PressureInterpolation());
}
/*}}}*/
int      Tetra::NumberofNodesVelocity(void){/*{{{*/
	return TetraRef::NumberofNodes(this->VelocityInterpolation());
}
/*}}}*/
int      Tetra::ObjectEnum(void){/*{{{*/

	return TetraEnum;

}/*}}}*/
int      Tetra::PressureInterpolation(void){/*{{{*/
	return TetraRef::PressureInterpolation(this->element_type);
}
/*}}}*/
void     Tetra::ReduceMatrices(ElementMatrix* Ke,ElementVector* pe){/*{{{*/

	if(pe){
		if(this->element_type==MINIcondensedEnum){
			int indices[3]={12,13,14};
			pe->StaticCondensation(Ke,3,&indices[0]);
		}
		else if(this->element_type==P1bubblecondensedEnum){
			int size   = nodes[4]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int offset = 0;
			for(int i=0;i<4;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int* indices=xNew<int>(size);
			for(int i=0;i<size;i++) indices[i] = offset+i;
			pe->StaticCondensation(Ke,size,indices);
			xDelete<int>(indices);
		}
	}

	if(Ke){
		if(this->element_type==MINIcondensedEnum){
			int indices[3]={12,13,14};
			Ke->StaticCondensation(3,&indices[0]);
		}
		else if(this->element_type==P1bubblecondensedEnum){
			int size   = nodes[4]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int offset = 0;
			for(int i=0;i<4;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int* indices=xNew<int>(size);
			for(int i=0;i<size;i++) indices[i] = offset+i;
			Ke->StaticCondensation(size,indices);
			xDelete<int>(indices);
		}
	}
}
/*}}}*/
void     Tetra::ResetFSBasalBoundaryCondition(void){/*{{{*/

	int numnodes = this->GetNumberOfNodes();

	int          approximation;
	IssmDouble*  vertexonbase= NULL;
	IssmDouble   slopex,slopey,groundedice;
	IssmDouble   xz_plane[6];

	/*For FS only: we want the CS to be tangential to the bedrock*/
	inputs->GetInputValue(&approximation,ApproximationEnum);
	if(!HasNodeOnBase() ||  approximation!=FSApproximationEnum) return;

	//printf("element number %i \n",this->id);
	/*Get inputs*/
	Input* slopex_input=inputs->GetInput(BedSlopeXEnum); _assert_(slopex_input);
	Input* slopey_input=inputs->GetInput(BedSlopeYEnum); _assert_(slopey_input);
	Input* groundedicelevelset_input=inputs->GetInput(MaskGroundediceLevelsetEnum); _assert_(groundedicelevelset_input);
	vertexonbase = xNew<IssmDouble>(numnodes);
	this->GetInputListOnNodesVelocity(&vertexonbase[0],MeshVertexonbaseEnum);

	/*Loop over basal nodes and update their CS*/
	GaussTetra* gauss = new GaussTetra();
	for(int i=0;i<this->NumberofNodesVelocity();i++){

		if(vertexonbase[i]==1){
			gauss->GaussNode(this->VelocityInterpolation(),i);

			slopex_input->GetInputValue(&slopex,gauss);
			slopey_input->GetInputValue(&slopey,gauss);
			groundedicelevelset_input->GetInputValue(&groundedice,gauss);

			/*New X axis          New Z axis*/
			xz_plane[0]=1.;       xz_plane[3]=-slopex;  
			xz_plane[1]=0.;       xz_plane[4]=-slopey;  
			xz_plane[2]=slopex;   xz_plane[5]=1.;          

			if(groundedice>0){
				if(this->nodes[i]->GetApproximation()==FSvelocityEnum){
					this->nodes[i]->DofInSSet(2); //vz 
				}
				else _error_("Flow equation approximation"<<EnumToStringx(this->nodes[i]->GetApproximation())<<" not supported yet");
			}
			else{
				if(this->nodes[i]->GetApproximation()==FSvelocityEnum){
					this->nodes[i]->DofInFSet(2); //vz
				}
				else _error_("Flow equation approximation"<<EnumToStringx(this->nodes[i]->GetApproximation())<<" not supported yet");
			}

			XZvectorsToCoordinateSystem(&this->nodes[i]->coord_system[0][0],&xz_plane[0]);
		}
	}

	/*cleanup*/
	xDelete<IssmDouble>(vertexonbase);
	delete gauss;
}
/*}}}*/
void     Tetra::ResetHooks(){/*{{{*/

	this->nodes=NULL;
	this->vertices=NULL;
	this->material=NULL;
	this->matpar=NULL;
	this->parameters=NULL;

	//deal with ElementHook mother class
	for(int i=0;i<this->numanalyses;i++) if(this->hnodes[i]) this->hnodes[i]->reset();
	this->hvertices->reset();
	this->hmaterial->reset();
	this->hmatpar->reset();
	if(this->hneighbors) this->hneighbors->reset();
}
/*}}}*/
void     Tetra::SetCurrentConfiguration(Elements* elementsin, Loads* loadsin, Nodes* nodesin, Materials* materialsin, Parameters* parametersin){/*{{{*/

	/*go into parameters and get the analysis_counter: */
	int analysis_counter;
	parametersin->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Get Element type*/
	this->element_type=this->element_type_list[analysis_counter];

	/*Pick up nodes*/
	if(this->hnodes[analysis_counter]) this->nodes=(Node**)this->hnodes[analysis_counter]->deliverp();
	else this->nodes=NULL;

}
/*}}}*/
Element* Tetra::SpawnBasalElement(void){/*{{{*/

	_assert_(HasFaceOnBase());

	int index1,index2,index3;
	this->FaceOnBaseIndices(&index1,&index2,&index3);
	return SpawnTria(index1,index2,index3);
}/*}}}*/
Element* Tetra::SpawnTopElement(void){/*{{{*/

	_assert_(HasFaceOnSurface());

	int index1,index2,index3;
	this->FaceOnSurfaceIndices(&index1,&index2,&index3);
	return SpawnTria(index1,index2,index3);
}/*}}}*/
Tria*    Tetra::SpawnTria(int index1,int index2,int index3){/*{{{*/

	int analysis_counter;

	/*go into parameters and get the analysis_counter: */
	this->parameters->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Create Tria*/
	Tria* tria=new Tria();
	tria->id=this->id;
	tria->inputs=(Inputs*)this->inputs->SpawnTriaInputs(index1,index2,index3);
	tria->parameters=this->parameters;
	tria->element_type=P1Enum; //Only P1 CG for now (TO BE CHANGED)
	this->SpawnTriaHook(xDynamicCast<ElementHook*>(tria),index1,index2,index3);

	/*Spawn material*/
	tria->material=(Material*)this->material->copy2(tria);

	/*recover nodes, material and matpar: */
	tria->nodes    = (Node**)tria->hnodes[analysis_counter]->deliverp();
	tria->vertices = (Vertex**)tria->hvertices->deliverp();
	tria->matpar   = (Matpar*)tria->hmatpar->delivers();

	/*Return new Tria*/
	return tria;
}
/*}}}*/
int      Tetra::TensorInterpolation(void){/*{{{*/
	return TetraRef::TensorInterpolation(this->element_type);
}
/*}}}*/
void     Tetra::Update(int index,IoModel* iomodel,int analysis_counter,int analysis_type,int finiteelement_type){ /*{{{*/

	/*Intermediaries*/
	int        i;
	int        tetra_vertex_ids[6];
	IssmDouble nodeinputs[6];
	IssmDouble yts;
	bool       dakota_analysis;
	bool       isFS;
	int        numnodes;
	int*       tetra_node_ids = NULL;

	/*Fetch parameters: */
	iomodel->FindConstant(&yts,"md.constants.yts");
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");
	iomodel->FindConstant(&isFS,"md.flowequation.isFS");

	/*Checks if debuging*/
	_assert_(iomodel->elements);

	/*Recover element type*/
	this->element_type_list[analysis_counter]=finiteelement_type;

	/*Recover vertices ids needed to initialize inputs*/
	for(i=0;i<4;i++) tetra_vertex_ids[i]=iomodel->elements[4*index+i]; //ids for vertices are in the elements array from Matlab

	/*Recover nodes ids needed to initialize the node hook.*/
	switch(finiteelement_type){
		case P1Enum:
			numnodes         = 4;
			tetra_node_ids   = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->nodecounter+iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->nodecounter+iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->nodecounter+iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->nodecounter+iomodel->elements[4*index+3];
			break;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			numnodes         = 5;
			tetra_node_ids   = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->nodecounter+iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->nodecounter+iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->nodecounter+iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->nodecounter+iomodel->elements[4*index+3];
			tetra_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+index+1;
			break;
		case P2Enum:
			numnodes        = 10;
			tetra_node_ids   = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->nodecounter+iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->nodecounter+iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->nodecounter+iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->nodecounter+iomodel->elements[4*index+3];
			tetra_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+0]+1;
			tetra_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+1]+1;
			tetra_node_ids[6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+2]+1;
			tetra_node_ids[7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+3]+1;
			tetra_node_ids[8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+4]+1;
			tetra_node_ids[9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+5]+1;
			break;
		case MINIEnum: case MINIcondensedEnum:
			numnodes         = 9;
			tetra_node_ids   = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->nodecounter+iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->nodecounter+iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->nodecounter+iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->nodecounter+iomodel->elements[4*index+3];
			tetra_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+index+1;

			tetra_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[4*index+0];
			tetra_node_ids[6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[4*index+1];
			tetra_node_ids[7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[4*index+2];
			tetra_node_ids[8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[4*index+3];
			break;
		case TaylorHoodEnum:
		case XTaylorHoodEnum:
			numnodes        = 14;
			tetra_node_ids  = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->nodecounter+iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->nodecounter+iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->nodecounter+iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->nodecounter+iomodel->elements[4*index+3];
			tetra_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+0]+1;
			tetra_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+1]+1;
			tetra_node_ids[6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+2]+1;
			tetra_node_ids[7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+3]+1;
			tetra_node_ids[8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+4]+1;
			tetra_node_ids[9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+5]+1;

			tetra_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[4*index+0];
			tetra_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[4*index+1];
			tetra_node_ids[12]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[4*index+2];
			tetra_node_ids[13]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[4*index+3];
			break;
		case LATaylorHoodEnum:
			numnodes        = 10;
			tetra_node_ids  = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->nodecounter+iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->nodecounter+iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->nodecounter+iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->nodecounter+iomodel->elements[4*index+3];
			tetra_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+0]+1;
			tetra_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+1]+1;
			tetra_node_ids[6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+2]+1;
			tetra_node_ids[7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+3]+1;
			tetra_node_ids[8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+4]+1;
			tetra_node_ids[9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+5]+1;
			break;
		default:
			_error_("Finite element "<<EnumToStringx(finiteelement_type)<<" not supported yet");
	}

	/*hooks: */
	this->SetHookNodes(tetra_node_ids,numnodes,analysis_counter); this->nodes=NULL;
	xDelete<int>(tetra_node_ids);

	/*Fill with IoModel*/
	this->InputUpdateFromIoModel(index,iomodel);
}
/*}}}*/
void     Tetra::ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss){/*{{{*/
	TetraRef::GetInputValue(pvalue,values,gauss,P1Enum);
}
/*}}}*/
int      Tetra::VelocityInterpolation(void){/*{{{*/
	return TetraRef::VelocityInterpolation(this->element_type);
}
/*}}}*/
void     Tetra::ViscousHeating(IssmDouble* pphi,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble phi;
	IssmDouble viscosity;
	IssmDouble epsilon[6];

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->StrainRateFS(&epsilon[0],xyz_list,(GaussTetra*)gauss,vx_input,vy_input,vz_input);
	this->material->ViscosityFS(&viscosity,3,xyz_list,(GaussTetra*)gauss,vx_input,vy_input,vz_input);
	GetPhi(&phi,&epsilon[0],viscosity);

	/*Assign output pointer*/
	*pphi = phi;
}
/*}}}*/
void     Tetra::ZeroLevelsetCoordinates(IssmDouble** pxyz_zero,IssmDouble* xyz_list,int levelsetenum){/*{{{*/
	/*Compute portion of the element that is grounded*/ 

	IssmDouble  levelset[NUMVERTICES];
	IssmDouble* xyz_zero = xNew<IssmDouble>(3*3);

	/*Recover parameters and values*/
	GetInputListOnVertices(&levelset[0],levelsetenum);

	if(levelset[0]==0. && levelset[1]==0. && levelset[2]==0.){ 
		xyz_zero[3*0+0]=xyz_list[0*3+0];
		xyz_zero[3*0+1]=xyz_list[0*3+1];
		xyz_zero[3*0+2]=xyz_list[0*3+2];

		/*New point 2*/
		xyz_zero[3*1+0]=xyz_list[1*3+0];
		xyz_zero[3*1+1]=xyz_list[1*3+1];
		xyz_zero[3*1+2]=xyz_list[1*3+2];

		/*New point 3*/
		xyz_zero[3*2+0]=xyz_list[2*3+0];
		xyz_zero[3*2+1]=xyz_list[2*3+1];
		xyz_zero[3*2+2]=xyz_list[2*3+2];
	}
	else if(levelset[0]==0. && levelset[1]==0. && levelset[3]==0.){ 
		xyz_zero[3*0+0]=xyz_list[0*3+0];
		xyz_zero[3*0+1]=xyz_list[0*3+1];
		xyz_zero[3*0+2]=xyz_list[0*3+2];

		/*New point 2*/
		xyz_zero[3*1+0]=xyz_list[1*3+0];
		xyz_zero[3*1+1]=xyz_list[1*3+1];
		xyz_zero[3*1+2]=xyz_list[1*3+2];

		/*New point 3*/
		xyz_zero[3*2+0]=xyz_list[3*3+0];
		xyz_zero[3*2+1]=xyz_list[3*3+1];
		xyz_zero[3*2+2]=xyz_list[3*3+2];
	}
	else if(levelset[1]==0. && levelset[2]==0. && levelset[3]==0.){ 
		xyz_zero[3*0+0]=xyz_list[1*3+0];
		xyz_zero[3*0+1]=xyz_list[1*3+1];
		xyz_zero[3*0+2]=xyz_list[1*3+2];

		/*New point 2*/
		xyz_zero[3*1+0]=xyz_list[2*3+0];
		xyz_zero[3*1+1]=xyz_list[2*3+1];
		xyz_zero[3*1+2]=xyz_list[2*3+2];

		/*New point 3*/
		xyz_zero[3*2+0]=xyz_list[3*3+0];
		xyz_zero[3*2+1]=xyz_list[3*3+1];
		xyz_zero[3*2+2]=xyz_list[3*3+2];
	}
	else if(levelset[2]==0. && levelset[0]==0. && levelset[3]==0.){ 
		xyz_zero[3*0+0]=xyz_list[2*3+0];
		xyz_zero[3*0+1]=xyz_list[2*3+1];
		xyz_zero[3*0+2]=xyz_list[2*3+2];

		/*New point 2*/
		xyz_zero[3*1+0]=xyz_list[0*3+0];
		xyz_zero[3*1+1]=xyz_list[0*3+1];
		xyz_zero[3*1+2]=xyz_list[0*3+2];

		/*New point 3*/
		xyz_zero[3*2+0]=xyz_list[3*3+0];
		xyz_zero[3*2+1]=xyz_list[3*3+1];
		xyz_zero[3*2+2]=xyz_list[3*3+2];
	}
	else _error_("Case not covered");

	/*Assign output pointer*/
	*pxyz_zero= xyz_zero;
}
/*}}}*/
