#include "./HydrologyDCEfficientAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
int  HydrologyDCEfficientAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/

void HydrologyDCEfficientAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int         hydrology_model;
	int         eplflip_lock;
	int         eplthickcomp;
	bool        isefficientlayer;
	/*retrieve some parameters: */
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want DC?*/
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");

	/*If not return*/
	if(!isefficientlayer) return;

	/*If yes, initialize a flip flop counter*/
	iomodel->FetchData(&eplflip_lock,"md.hydrology.eplflip_lock");
	parameters->AddObject(new IntParam(HydrologydcEplflipLockEnum,eplflip_lock));

	iomodel->FetchData(&eplthickcomp,"md.hydrology.epl_thick_comp");
	parameters->AddObject(new IntParam(HydrologydcEplThickCompEnum,eplthickcomp));
}/*}}}*/

void HydrologyDCEfficientAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	bool   isefficientlayer;
	int    hydrology_model;

	/*Now, do we really want DC?*/
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
	if(!isefficientlayer) return;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}
	iomodel->FetchDataToInput(elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.epl_head",EplHeadEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.sediment_head",SedimentHeadEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.epl_thickness",HydrologydcEplThicknessEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.basal_moulin_input",HydrologydcBasalMoulinInputEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
}/*}}}*/

void HydrologyDCEfficientAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	/*Now, do we really want DC?*/
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	bool isefficientlayer;
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
	if(!isefficientlayer) return;

	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	}
	::CreateNodes(nodes,iomodel,HydrologyDCEfficientAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/

void HydrologyDCEfficientAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Do we really want DC?*/
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	bool isefficientlayer;
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
	if(!isefficientlayer) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spcepl_head",HydrologyDCEfficientAnalysisEnum,P1Enum);
}/*}}}*/

void HydrologyDCEfficientAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Do we really want DC?*/
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	/*Do we want an efficient layer*/
	bool isefficientlayer;
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
	if(!isefficientlayer) return;

	/*Fetch parameters: */
	if(iomodel->domaintype==Domain3DEnum){
		iomodel->FetchData(1,"md.mesh.vertexonbase");
	}

	//Add moulin inputs as loads
	CreateSingleNodeToElementConnectivity(iomodel);
	for(int i=0;i<iomodel->numberofvertices;i++){
		if (iomodel->domaintype!=Domain3DEnum){
			/*keep only this partition's nodes:*/
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(iomodel->loadcounter+i+1,i,iomodel,HydrologyDCEfficientAnalysisEnum));
			}
		}
		else if(reCast<int>(iomodel->Data("md.mesh.vertexonbase")[i])){
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(iomodel->loadcounter+i+1,i,iomodel,HydrologyDCEfficientAnalysisEnum));
			}	
		}
	}
	iomodel->DeleteData(1,"md.mesh.vertexonbase");
}/*}}}*/

void HydrologyDCEfficientAnalysis::InitZigZagCounter(FemModel* femmodel){/*{{{*/

	int*   eplzigzag_counter =NULL;
	eplzigzag_counter=xNewZeroInit<int>(femmodel->nodes->Size());
	femmodel->parameters->AddObject(new IntVecParam(EplZigZagCounterEnum,eplzigzag_counter,femmodel->nodes->Size()));
	xDelete<int>(eplzigzag_counter);
}/*}}}*/

void HydrologyDCEfficientAnalysis::ResetCounter(FemModel* femmodel){/*{{{*/

	int*     eplzigzag_counter=NULL;
	femmodel->parameters->FindParam(&eplzigzag_counter,NULL,EplZigZagCounterEnum);
	for(int i=0;i<femmodel->nodes->Size();i++){
		eplzigzag_counter[i]=0;
	}
	femmodel->parameters->SetParam(eplzigzag_counter,femmodel->nodes->Size(),EplZigZagCounterEnum);
	xDelete<int>(eplzigzag_counter);
}/*}}}*/

/*Finite Element Analysis*/
void HydrologyDCEfficientAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/

ElementVector* HydrologyDCEfficientAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/

ElementMatrix* HydrologyDCEfficientAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/

ElementMatrix* HydrologyDCEfficientAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	bool     active_element;
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	Input* active_element_input = basalelement->GetInput(HydrologydcMaskEplactiveEltEnum); _assert_(active_element_input);
	active_element_input->GetInputValue(&active_element);

	/*Check that all nodes are active, else return empty matrix*/
	if(!active_element) {
	if(domaintype!=Domain2DhorizontalEnum){
			basalelement->DeleteMaterials(); 
			delete basalelement;
		}
		return NULL;
	}

	/* Intermediaries */
	IssmDouble  D_scalar,Jdet,dt;
	IssmDouble  transfer;
	IssmDouble  epl_transmitivity;
	IssmDouble  epl_storing;
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke     = basalelement->NewElementMatrix();
	IssmDouble*    B      = xNew<IssmDouble>(2*numnodes);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble     D[2][2]={0.};

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);

	Input* epl_thick_input = basalelement->GetInput(HydrologydcEplThicknessEnum); _assert_(epl_thick_input);
	Input* epl_head_input	= basalelement->GetInput(EplHeadEnum);  _assert_(epl_head_input);
	Input* base_input			= basalelement->GetInput(BaseEnum); _assert_(base_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss					= basalelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss           ->GaussPoint(ig);
		basalelement    ->JacobianDeterminant(&Jdet,xyz_list,gauss);
		
		epl_transmitivity = EplTransmitivity(basalelement,gauss,epl_thick_input,epl_head_input,base_input);
		epl_storing				= EplStoring(basalelement,gauss,epl_thick_input,epl_head_input,base_input);

		/*Diffusivity*/
		D_scalar=epl_transmitivity*gauss->weight*Jdet;
		if(dt!=0.) D_scalar=D_scalar*dt;
		D[0][0]=D_scalar;
		D[1][1]=D_scalar;
		GetB(B,basalelement,xyz_list,gauss); 
		TripleMultiply(B,2,numnodes,1,
					&D[0][0],2,2,0,
					B,2,numnodes,0,
					&Ke->values[0],1);

		/*Transient*/
		if(dt!=0.){
			basalelement->NodalFunctions(&basis[0],gauss);
			D_scalar=epl_storing*gauss->weight*Jdet;
			TripleMultiply(basis,numnodes,1,0,
						&D_scalar,1,1,0,
						basis,1,numnodes,0,
						&Ke->values[0],1);
			
			
			/*Transfer EPL part*/
			transfer=GetHydrologyKMatrixTransfer(basalelement);
			D_scalar=dt*transfer*gauss->weight*Jdet;
			TripleMultiply(basis,numnodes,1,0,
										 &D_scalar,1,1,0,
										 basis,1,numnodes,0,
										 &Ke->values[0],1);
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(B);
	delete gauss;
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/

ElementVector* HydrologyDCEfficientAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
	bool     active_element;
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	Input* active_element_input = basalelement->GetInput(HydrologydcMaskEplactiveEltEnum); _assert_(active_element_input);
	active_element_input->GetInputValue(&active_element);

	/*Check that all nodes are active, else return empty matrix*/
	if(!active_element) {
		if(domaintype!=Domain2DhorizontalEnum){
			basalelement->DeleteMaterials(); 
			delete basalelement;
		}
		return NULL;
	}
	/*Intermediaries */
	IssmDouble dt,scalar,water_head;
	IssmDouble water_load,transfer;
	IssmDouble epl_storing;
	IssmDouble Jdet;
	IssmDouble residual,connectivity;

	IssmDouble *xyz_list     = NULL;
	Input*      old_wh_input = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = basalelement->GetNumberOfNodes();
	int numvertices = basalelement->GetNumberOfVertices();

	/*Initialize Element vector*/
	ElementVector* pe    = basalelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);	

	Input* epl_thick_input = basalelement->GetInput(HydrologydcEplThicknessEnum); _assert_(epl_thick_input);
	Input* sed_head_input  = basalelement->GetInput(SedimentHeadEnum); _assert_(sed_head_input);
	Input* epl_head_input	 = basalelement->GetInput(EplHeadEnum); _assert_(epl_head_input);
	Input* water_input		 = basalelement->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(water_input);
	Input* residual_input  = basalelement->GetInput(SedimentHeadResidualEnum); _assert_(residual_input);
	Input* base_input			 = basalelement->GetInput(BaseEnum); _assert_(base_input);

	if(dt!= 0.){
		old_wh_input = basalelement->GetInput(EplHeadOldEnum);            _assert_(old_wh_input);
	}
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss           = basalelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);
		basalelement ->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement ->NodalFunctions(basis,gauss);
		epl_storing	= EplStoring(basalelement,gauss,epl_thick_input,epl_head_input,base_input);
		/*Loading term*/
		water_input->GetInputValue(&water_load,gauss);
		scalar = Jdet*gauss->weight*(water_load);
		if(dt!=0.) scalar = scalar*dt;
		for(int i=0;i<numnodes;i++)pe->values[i]+=scalar*basis[i];
		
		/*Transient and transfer terms*/
		if(dt!=0.){
			old_wh_input->GetInputValue(&water_head,gauss);
			/*Dealing with the epl part of the transfer term*/
			transfer=GetHydrologyPVectorTransfer(basalelement,gauss,sed_head_input);
			scalar = Jdet*gauss->weight*((water_head*epl_storing)+(dt*transfer));
			for(int i=0;i<numnodes;i++)pe->values[i]+=scalar*basis[i];
		}
	}
	delete gauss;

	/*	Add residual if necessary*/
	gauss=basalelement->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);
		connectivity = IssmDouble(basalelement->VertexConnectivity(iv));
		residual_input->GetInputValue(&residual,gauss);
		pe->values[iv]+=residual/connectivity;
	}
	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/

void HydrologyDCEfficientAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,EplHeadEnum);
}/*}}}*/

void HydrologyDCEfficientAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/

void HydrologyDCEfficientAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	/*Intermediaries*/
	bool     active_element;
	int      domaintype;
	Element* basalelement=NULL;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	/*Intermediary*/
	int* doflist = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	IssmDouble* eplHeads    = xNew<IssmDouble>(numnodes);
	basalelement->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);

	Input* active_element_input=basalelement->GetInput(HydrologydcMaskEplactiveEltEnum); _assert_(active_element_input);
	active_element_input->GetInputValue(&active_element);

	/*Use the dof list to index into the solution vector: */
	/*If the EPL is not active we revert to the bedrock elevation*/
	if(active_element){
		for(int i=0;i<numnodes;i++){
			eplHeads[i]=solution[doflist[i]];
			if(xIsNan<IssmDouble>(eplHeads[i])) _error_("NaN found in solution vector");
			if(xIsInf<IssmDouble>(eplHeads[i])) _error_("Inf found in solution vector");
		}
	}
	else{
		basalelement->GetInputListOnVertices(&eplHeads[0],BaseEnum);
		for(int i=0;i<numnodes;i++){
			if(xIsNan<IssmDouble>(eplHeads[i])) _error_("NaN found in solution vector");
			if(xIsInf<IssmDouble>(eplHeads[i])) _error_("Inf found in solution vector");
		}
	}
	/*Add input to the element: */
	element->AddBasalInput(EplHeadEnum,eplHeads,P1Enum);
	/*Free ressources:*/
	xDelete<IssmDouble>(eplHeads);
	xDelete<int>(doflist);
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
} /*}}}*/

void HydrologyDCEfficientAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/

/*Intermediaries*/
IssmDouble HydrologyDCEfficientAnalysis::GetHydrologyKMatrixTransfer(Element* element){/*{{{*/

	int transfermethod;
	IssmDouble leakage,transfer;

	element->FindParam(&transfermethod,HydrologydcTransferFlagEnum);
	/*Switch between the different transfer methods cases*/
	switch(transfermethod){
	case 0:
		/*Just keepping the transfer to zero*/
		transfer=0.0;
		break;
	case 1:
		element->FindParam(&leakage,HydrologydcLeakageFactorEnum);
		transfer=leakage;
		break;
	default:
		_error_("no case higher than 1 for the Transfer method");
	}
	return transfer;
}/*}}}*/

IssmDouble HydrologyDCEfficientAnalysis::GetHydrologyPVectorTransfer(Element* element, Gauss* gauss, Input* sed_head_input){/*{{{*/

	int transfermethod;
	IssmDouble sediment_head;
	IssmDouble leakage,transfer;

	element->FindParam(&transfermethod,HydrologydcTransferFlagEnum);
	/*Switch between the different transfer methods cases*/
	switch(transfermethod){
	case 0:
		/*Just keepping the transfer to zero*/
		transfer=0.0;
		break;
	case 1:
		_assert_(sed_head_input);
		sed_head_input->GetInputValue(&sediment_head,gauss);
		element->FindParam(&leakage,HydrologydcLeakageFactorEnum);
		transfer=sediment_head*leakage;
		break;
	default:
		_error_("no case higher than 1 for the Transfer method");
	}

	return transfer;
}/*}}}*/

void HydrologyDCEfficientAnalysis::ComputeEPLThickness(FemModel* femmodel){/*{{{*/

	bool        active_element;
	int         iseplthickcomp;
	int         domaintype;
	IssmDouble  dt,A,B;
	IssmDouble  EPLgrad2;
	IssmDouble  EPL_N;
	
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&iseplthickcomp,HydrologydcEplThickCompEnum);
	if(iseplthickcomp==0) return;

	for(int j=0;j<femmodel->elements->Size();j++){
		
		Element* element=(Element*)femmodel->elements->GetObjectByOffset(j);
		
		switch(domaintype){
		case Domain2DhorizontalEnum:
			if(!element->IsOnBase()) return;	
			B = element->GetMaterialParameter(MaterialsRheologyBbarEnum);
			break;
		case Domain3DEnum:
			B = element->GetMaterialParameter(MaterialsRheologyBEnum);
			break;
		default:
		_error_("not Implemented Yet");
		}
		
		int         numnodes      = element->GetNumberOfNodes();
		IssmDouble* thickness     = xNew<IssmDouble>(numnodes);
		IssmDouble* eplhead       = xNew<IssmDouble>(numnodes);
		IssmDouble* epl_slopeX    = xNew<IssmDouble>(numnodes);
		IssmDouble* epl_slopeY    = xNew<IssmDouble>(numnodes);
		IssmDouble* old_thickness = xNew<IssmDouble>(numnodes);
		IssmDouble* ice_thickness = xNew<IssmDouble>(numnodes);
		IssmDouble* bed           = xNew<IssmDouble>(numnodes);

		Input* 	active_element_input=element->GetInput(HydrologydcMaskEplactiveEltEnum); _assert_(active_element_input);		
		active_element_input->GetInputValue(&active_element);
		element->FindParam(&dt,TimesteppingTimeStepEnum);
	
		/*For now, assuming just one way to compute EPL thickness*/
		IssmDouble gravity          = element->GetMaterialParameter(ConstantsGEnum);
		IssmDouble rho_water        = element->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
		IssmDouble rho_ice          = element->GetMaterialParameter(MaterialsRhoIceEnum);
		IssmDouble n                =	element->GetMaterialParameter(MaterialsRheologyNEnum);
		IssmDouble latentheat       = element->GetMaterialParameter(MaterialsLatentheatEnum);
		IssmDouble epl_conductivity = element->GetMaterialParameter(HydrologydcEplConductivityEnum);
		IssmDouble init_thick       =	element->GetMaterialParameter(HydrologydcEplInitialThicknessEnum);
		IssmDouble max_thick        =	element->GetMaterialParameter(HydrologydcEplMaxThicknessEnum);
		
		A=pow(B,-n);
		
		element->GetInputListOnVertices(&eplhead[0],EplHeadEnum);
		element->GetInputListOnVertices(&epl_slopeX[0],EplHeadSlopeXEnum); 
		element->GetInputListOnVertices(&epl_slopeY[0],EplHeadSlopeYEnum);
		element->GetInputListOnVertices(&old_thickness[0],HydrologydcEplThicknessOldEnum);
		element->GetInputListOnVertices(&ice_thickness[0],ThicknessEnum);
		element->GetInputListOnVertices(&bed[0],BaseEnum);
		
		if(!active_element){
			/*Keeping thickness to initial value if EPL is not active*/
			for(int i=0;i<numnodes;i++){
				thickness[i]=init_thick;
			}
		}
		else{
			for(int i=0;i<numnodes;i++){
				/*Compute first the effective pressure in the EPL*/
				EPL_N=gravity*((rho_ice*ice_thickness[i])-(rho_water*(eplhead[i]-bed[i])));
				if(EPL_N<0.0)EPL_N=0.0;
				/*Get then the square of the gradient of EPL heads*/
				EPLgrad2 = (epl_slopeX[i]*epl_slopeX[i])+(epl_slopeY[i]*epl_slopeY[i]);
				/*And proceed to the real thing*/
				thickness[i] = old_thickness[i]/(1.0
																				 -((rho_water*gravity*epl_conductivity*EPLgrad2*dt)/(rho_ice*latentheat))
																				 +((2.0*A*dt*pow(EPL_N,n))/(pow(n,n))));
				/*Take care of otherthikening*/
				if(thickness[i]>max_thick){
					thickness[i] = max_thick;
				}
			}
		}
		element->AddInput(HydrologydcEplThicknessEnum,thickness,element->GetElementType());
		xDelete<IssmDouble>(thickness);
		xDelete<IssmDouble>(eplhead);
		xDelete<IssmDouble>(epl_slopeX);
		xDelete<IssmDouble>(epl_slopeY);
		xDelete<IssmDouble>(old_thickness);
		xDelete<IssmDouble>(ice_thickness);
		xDelete<IssmDouble>(bed);
	}
}/*}}}*/

void HydrologyDCEfficientAnalysis::GetB(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3] where Bi is of size 3*NDOF2. 
	 * For node i, Bi can be expressed in the actual coordinate system
	 * by: 
	 *       Bi=[ dN/dx ]
	 *          [ dN/dy ]
	 * where N is the finiteelement function for node i.
	 *
	 * We assume B has been allocated already, of size: 3x(NDOF2*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* dbasis=xNew<IssmDouble>(2*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B: */
	for(int i=0;i<numnodes;i++){
		B[numnodes*0+i] = dbasis[0*numnodes+i];
		B[numnodes*1+i] = dbasis[1*numnodes+i];
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);
}/*}}}*/

void  HydrologyDCEfficientAnalysis::HydrologyEPLGetMask(Vector<IssmDouble>* vec_mask, Vector<IssmDouble>* recurence, int* eplzigzag_counter, Element* element){

	bool        active_element;
	int         domaintype;
	IssmDouble  h_max;
	IssmDouble  sedheadmin;
	Element*    basalelement=NULL;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	/*Intermediaries*/
	int         numnodes      =basalelement->GetNumberOfNodes();
	IssmDouble* epl_thickness =xNew<IssmDouble>(numnodes);
	IssmDouble* old_active    =xNew<IssmDouble>(numnodes);
	IssmDouble* sedhead       =xNew<IssmDouble>(numnodes);
	IssmDouble* eplhead       =xNew<IssmDouble>(numnodes);
	IssmDouble* residual      =xNew<IssmDouble>(numnodes);
	IssmDouble* base          =xNew<IssmDouble>(numnodes);
	
	IssmDouble init_thick    =basalelement->GetMaterialParameter(HydrologydcEplInitialThicknessEnum);
	IssmDouble colapse_thick =basalelement->GetMaterialParameter(HydrologydcEplColapseThicknessEnum);

	Input* active_element_input=basalelement->GetInput(HydrologydcMaskEplactiveEltEnum); _assert_(active_element_input);
	active_element_input->GetInputValue(&active_element);

	basalelement-> GetInputListOnVertices(&old_active[0],HydrologydcMaskEplactiveNodeEnum);	
	basalelement-> GetInputListOnVertices(&epl_thickness[0],HydrologydcEplThicknessEnum);	
	basalelement-> GetInputListOnVertices(&sedhead[0],SedimentHeadEnum);
	basalelement-> GetInputListOnVertices(&eplhead[0],EplHeadEnum);
	basalelement-> GetInputListOnVertices(&residual[0],SedimentHeadResidualEnum);
	basalelement-> GetInputListOnVertices(&base[0],BaseEnum);

	/*Get minimum sediment head of the element*/
	sedheadmin=sedhead[0];
	for(int i=1;i<numnodes;i++) if(sedhead[i]<=sedheadmin)sedheadmin=sedhead[i];
	for(int i=0;i<numnodes;i++){
		GetHydrologyDCInefficientHmax(&h_max,basalelement,basalelement->nodes[i]);
		/*If mask was already one, keep one or colapse*/
		if(old_active[i]>0.){
			vec_mask->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
			/* If epl thickness gets under colapse thickness, close the layer */
			if(epl_thickness[i]<colapse_thick){
				vec_mask->SetValue(basalelement->nodes[i]->Sid(),0.,INS_VAL);
				recurence->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
			}
			/* //If epl head gets under base elevation, close the layer */
			/* else if(eplhead[i]<(base[i]-1.0e-8)){ */
			/* 	vec_mask->SetValue(basalelement->nodes[i]->Sid(),0.,INS_VAL); */
			/* 	recurence->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL); */
			/* } */
		}
		/*If node is now closed bring its thickness back to initial*/
		if (old_active[i]==0.){
			epl_thickness[i]=init_thick;
			/*Activate if we have a residual from sediment*/
			if(residual[i]>0.){
				vec_mask->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
				if(old_active[i]==0.){
					recurence->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
				}
			}
		}
		/*Increase of the efficient system is needed if the epl head reach the maximum value (sediment max value for now)*/
		if(eplhead[i]>=h_max && active_element){
			for(int j=0;j<numnodes;j++){
				/*Increase of the domain is on the downstream node in term of sediment head*/
				if((sedhead[j] == sedheadmin) && (i!=j)){
					vec_mask->SetValue(basalelement->nodes[j]->Sid(),1.,INS_VAL);
					if(old_active[j]==0.){
						recurence->SetValue(basalelement->nodes[j]->Sid(),1.,INS_VAL);
					}
				}
			}
		}
	}
	basalelement->AddInput(HydrologydcEplThicknessEnum,epl_thickness,basalelement->GetElementType());

	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	xDelete<IssmDouble>(epl_thickness);
	xDelete<IssmDouble>(old_active);
	xDelete<IssmDouble>(sedhead);
	xDelete<IssmDouble>(eplhead);
	xDelete<IssmDouble>(residual);
	xDelete<IssmDouble>(base);
}
/*}}}*/
IssmDouble HydrologyDCEfficientAnalysis::EplStoring(Element* element,Gauss* gauss, Input* epl_thick_input, Input* epl_head_input, Input* base_input){/*{{{*/
	IssmDouble epl_storing;
	IssmDouble water_sheet,storing;
	IssmDouble epl_thickness,prestep_head,base_elev;
	IssmDouble rho_freshwater        = element->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	IssmDouble g                     = element->GetMaterialParameter(ConstantsGEnum);
	IssmDouble epl_porosity					 = element->GetMaterialParameter(HydrologydcEplPorosityEnum);
	IssmDouble epl_compressibility	 = element->GetMaterialParameter(HydrologydcEplCompressibilityEnum);
	IssmDouble water_compressibility = element->GetMaterialParameter(HydrologydcWaterCompressibilityEnum);

	epl_thick_input->GetInputValue(&epl_thickness,gauss);
	epl_head_input->GetInputValue(&prestep_head,gauss);
	base_input->GetInputValue(&base_elev,gauss);
	water_sheet=max(0.0,(prestep_head-base_elev));

	storing=rho_freshwater*g*epl_porosity*epl_thickness*(water_compressibility+(epl_compressibility/epl_porosity));

	/* //porosity for unconfined region */
	/* if (water_sheet<=0.9*epl_thickness){ */
	/* 	epl_storing=epl_porosity; */
	/* } */
	/* //continuity ramp */
	/* else if((water_sheet<epl_thickness) && (water_sheet>0.9*epl_thickness)){ */
	/* 	epl_storing=(epl_thickness-water_sheet)*(epl_porosity-storing)/(0.1*epl_thickness)+storing; */
	/* } */
	/* //storing coefficient for confined */
	/* else{ */
	/* 	epl_storing=storing; */
	/* } */
 	/* return epl_storing; */
	return storing;
}/*}}}*/

IssmDouble HydrologyDCEfficientAnalysis::EplTransmitivity(Element* element,Gauss* gauss, Input* epl_thick_input, Input* epl_head_input, Input* base_input){/*{{{*/
	IssmDouble epl_transmitivity;
	IssmDouble water_sheet;
	IssmDouble epl_thickness,base_elev,prestep_head;
	IssmDouble epl_conductivity      = element->GetMaterialParameter(HydrologydcEplConductivityEnum);
	epl_thick_input->GetInputValue(&epl_thickness,gauss);
	epl_head_input->GetInputValue(&prestep_head,gauss);
	base_input->GetInputValue(&base_elev,gauss);

	water_sheet=max(0.0,(prestep_head-base_elev));
	
	epl_transmitivity=epl_conductivity*epl_thickness;
	//epl_transmitivity=max(1.0e-6,(epl_conductivity*min(water_sheet,epl_thickness)));
	return epl_transmitivity;
}/*}}}*/

void HydrologyDCEfficientAnalysis::HydrologyEPLGetActive(Vector<IssmDouble>* active_vec, Element* element){/*{{{*/
	/*Constants*/
	int      domaintype;
	Element*   basalelement=NULL;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	
	const int   numnodes = basalelement->GetNumberOfNodes();
	IssmDouble  flag     = 0.;
	IssmDouble* active   = xNew<IssmDouble>(numnodes);
	bool active_element;

	/*Pass the activity mask from elements to nodes*/
	basalelement->GetInputListOnVertices(&active[0],HydrologydcMaskEplactiveNodeEnum);
	Input* 	active_element_input=basalelement->GetInput(HydrologydcMaskEplactiveEltEnum); _assert_(active_element_input);		
	active_element_input->GetInputValue(&active_element);
	
	for(int i=0;i<numnodes;i++) flag+=active[i];

	/*If any node is active all the node in the element are active*/
	if(flag>0.){
		for(int i=0;i<numnodes;i++){
			active_vec->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
		}
	}
	/*If the element is active all its nodes are active*/
	else if(active_element){
		for(int i=0;i<numnodes;i++){
			active_vec->SetValue(basalelement->nodes[i]->Sid(),1.,INS_VAL);
		}		
	}
	else{
		/*Do not do anything: at least one node is active for this element but this element is not solved for*/
	}
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	xDelete<IssmDouble>(active);
}
/*}}}*/

void HydrologyDCEfficientAnalysis::GetHydrologyDCInefficientHmax(IssmDouble* ph_max,Element* element, Node* innode){/*{{{*/
	
	int        hmax_flag;
	IssmDouble h_max;
	IssmDouble rho_ice,rho_water;
	IssmDouble thickness,bed;
	/*Get the flag to the limitation method*/
	element->FindParam(&hmax_flag,HydrologydcSedimentlimitFlagEnum);
	
	/*Switch between the different cases*/
	switch(hmax_flag){
	case 0:
		h_max=1.0e+10;
		break;
	case 1:
		element->FindParam(&h_max,HydrologydcSedimentlimitEnum);
		break;
	case 2:
		/*Compute max*/
		rho_water = element->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
		rho_ice   = element->GetMaterialParameter(MaterialsRhoIceEnum);
		element-> GetInputValue(&thickness,innode,ThicknessEnum);
		element-> GetInputValue(&bed,innode,BaseEnum);
		h_max=((rho_ice*thickness)/rho_water)+bed;
		break;
	case 3:
		_error_("Using normal stress  not supported yet");
		break;
	default:
		_error_("no case higher than 3 for SedimentlimitFlag");
	}
	/*Assign output pointer*/
	*ph_max=h_max;
}
/*}}}*/
