#include "./HydrologyDCInefficientAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../classes/Node.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
int  HydrologyDCInefficientAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/

void HydrologyDCInefficientAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int         hydrology_model;
	int         sedimentlimit_flag;
	int         transfer_flag;
	int         penalty_lock;
	int         hydro_maxiter;
	bool        isefficientlayer;
	IssmDouble  penalty_factor;
	IssmDouble  rel_tol;
	IssmDouble  leakagefactor;
	IssmDouble  sedimentlimit;

	/*retrieve some parameters: */
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want DC?*/
	if(hydrology_model!=HydrologydcEnum) return;

	iomodel->FetchData(&sedimentlimit_flag, "md.hydrology.sedimentlimit_flag" );
	iomodel->FetchData(&transfer_flag,      "md.hydrology.transfer_flag" );
	iomodel->FetchData(&penalty_factor,     "md.hydrology.penalty_factor" );
	iomodel->FetchData(&isefficientlayer,   "md.hydrology.isefficientlayer");
	iomodel->FetchData(&hydro_maxiter,      "md.hydrology.max_iter" );
	iomodel->FetchData(&penalty_lock,       "md.hydrology.penalty_lock" );
	iomodel->FetchData(&rel_tol,            "md.hydrology.rel_tol" );

	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));
	parameters->AddObject(new IntParam(HydrologydcSedimentlimitFlagEnum,sedimentlimit_flag));
	parameters->AddObject(new IntParam(HydrologydcTransferFlagEnum,transfer_flag));
	parameters->AddObject(new IntParam(HydrologydcPenaltyLockEnum,penalty_lock));
	parameters->AddObject(new IntParam(HydrologydcMaxIterEnum,hydro_maxiter));
	parameters->AddObject(new BoolParam(HydrologydcIsefficientlayerEnum,isefficientlayer));
	parameters->AddObject(new DoubleParam(HydrologydcPenaltyFactorEnum,penalty_factor));
	parameters->AddObject(new DoubleParam(HydrologydcRelTolEnum,rel_tol));
	if(transfer_flag==1){
		iomodel->FetchData(&leakagefactor,"md.hydrology.leakage_factor");
		parameters->AddObject(new DoubleParam(HydrologydcLeakageFactorEnum,leakagefactor));
	}
	if(sedimentlimit_flag==1){
		iomodel->FetchData(&sedimentlimit,"md.hydrology.sedimentlimit");
		parameters->AddObject(new DoubleParam(HydrologydcSedimentlimitEnum,sedimentlimit));
	}
}/*}}}*/

void HydrologyDCInefficientAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	bool   isefficientlayer;
	bool   element_active;
	int    hydrology_model;
	
	/*Fetch data needed: */
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want DC?*/
	if(hydrology_model!=HydrologydcEnum) return;

	/*Fetch data needed: */
	iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");

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
	iomodel->FetchDataToInput(elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.basal_moulin_input",HydrologydcBasalMoulinInputEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.sediment_head",SedimentHeadEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.sediment_transmitivity",HydrologydcSedimentTransmitivityEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}

	if(isefficientlayer){
		iomodel->FetchDataToInput(elements,"md.hydrology.mask_eplactive_node",HydrologydcMaskEplactiveNodeEnum);
	}
}/*}}}*/

void HydrologyDCInefficientAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want DC?*/
	if(hydrology_model!=HydrologydcEnum) return;

	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	}
	::CreateNodes(nodes,iomodel,HydrologyDCInefficientAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/

void HydrologyDCInefficientAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spcsediment_head",HydrologyDCInefficientAnalysisEnum,P1Enum);
}/*}}}*/

void HydrologyDCInefficientAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	if(hydrology_model!=HydrologydcEnum) return;

	if(iomodel->domaintype==Domain3DEnum){
		iomodel->FetchData(1,"md.mesh.vertexonbase");
	}
	//create penalties for nodes: no node can have water above the max
	CreateSingleNodeToElementConnectivity(iomodel);
	for(int i=0;i<iomodel->numberofvertices;i++){
		if (iomodel->domaintype!=Domain3DEnum){
			/*keep only this partition's nodes:*/
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Pengrid(iomodel->loadcounter+i+1,i,iomodel,HydrologyDCInefficientAnalysisEnum));
				loads->AddObject(new Moulin(iomodel->loadcounter+i+1,i,iomodel,HydrologyDCInefficientAnalysisEnum));
			}
		}
		else if(reCast<int>(iomodel->Data("md.mesh.vertexonbase")[i])){
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Pengrid(iomodel->loadcounter+i+1,i,iomodel,HydrologyDCInefficientAnalysisEnum));
				loads->AddObject(new Moulin(iomodel->loadcounter+i+1,i,iomodel,HydrologyDCInefficientAnalysisEnum));
			}	
		}
	}
	iomodel->DeleteData(1,"md.mesh.vertexonbase");
}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyDCInefficientAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/

ElementVector* HydrologyDCInefficientAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/

ElementMatrix* HydrologyDCInefficientAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/

ElementMatrix* HydrologyDCInefficientAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
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

	/*Intermediaries */
	bool        active_element,isefficientlayer;
	IssmDouble  D_scalar,Jdet,dt;
	IssmDouble  sediment_transmitivity;
	IssmDouble  transfer,sediment_storing;
	IssmDouble *xyz_list  = NULL;

	/*Define transfer related variables*/
	Input* active_element_input =NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke     = basalelement->NewElementMatrix();
	IssmDouble*    B      = xNew<IssmDouble>(2*numnodes);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble     D[2][2]= {0.};

	/*Retrieve all inputs and parameters*/
	basalelement ->GetVerticesCoordinates(&xyz_list);
	basalelement ->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement ->FindParam(&isefficientlayer,HydrologydcIsefficientlayerEnum);
	Input* SedTrans_input = basalelement->GetInput(HydrologydcSedimentTransmitivityEnum); _assert_(SedTrans_input);
	Input* sed_head_input = basalelement->GetInput(SedimentHeadEnum);
	Input* base_input     = basalelement->GetInput(BaseEnum);

	/*Transfer related Inputs*/
	if(isefficientlayer){
		active_element_input = basalelement->GetInput(HydrologydcMaskEplactiveEltEnum); _assert_(active_element_input);
	}
	
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	for(int ig=gauss -> begin();ig<gauss->end();ig++){
		gauss          -> GaussPoint(ig);
		basalelement   -> JacobianDeterminant(&Jdet,xyz_list,gauss);

		sediment_transmitivity = SedimentTransmitivity(basalelement,gauss,sed_head_input,base_input,SedTrans_input);
		sediment_storing       = SedimentStoring(basalelement,gauss,sed_head_input,base_input);

		/*Diffusivity*/
		D_scalar=sediment_transmitivity*gauss->weight*Jdet;
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
			D_scalar=sediment_storing*gauss->weight*Jdet;
			TripleMultiply(basis,numnodes,1,0,
										 &D_scalar,1,1,0,
										 basis,1,numnodes,0,
										 &Ke->values[0],1);
			
			/*Transfer EPL part*/
			if(isefficientlayer){
				active_element_input->GetInputValue(&active_element);
				if(active_element){
					transfer=GetHydrologyKMatrixTransfer(basalelement);
					basalelement->NodalFunctions(&basis[0],gauss);
					D_scalar=dt*transfer*gauss->weight*Jdet;
					TripleMultiply(basis,numnodes,1,0,
												 &D_scalar,1,1,0,
												 basis,1,numnodes,0,
												 &Ke->values[0],1);
				}
			}
		}
	}
	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/

ElementVector* HydrologyDCInefficientAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
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

	/*Intermediaries */
	bool       active_element,isefficientlayer;
	IssmDouble dt,scalar,sediment_storing;
	IssmDouble water_head;
	IssmDouble water_load,transfer;
	IssmDouble Jdet;

	IssmDouble *xyz_list             = NULL;
	Input*      active_element_input = NULL;
	Input*      old_wh_input = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe    = basalelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement->FindParam(&isefficientlayer,HydrologydcIsefficientlayerEnum);

	Input* sed_head_input = basalelement->GetInput(SedimentHeadEnum);
	Input* epl_head_input = basalelement->GetInput(EplHeadEnum);
	Input* base_input		  = basalelement->GetInput(BaseEnum);
	Input* water_input	  = basalelement->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(water_input);
	if(dt!= 0.){
		old_wh_input = basalelement->GetInput(SedimentHeadOldEnum);                  _assert_(old_wh_input);
	}
	/*Transfer related Inputs*/
	if(isefficientlayer){
		active_element_input = basalelement->GetInput(HydrologydcMaskEplactiveEltEnum); _assert_(active_element_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);
		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		/*Loading term*/
		if(!isefficientlayer){
			water_input->GetInputValue(&water_load,gauss);
			scalar = Jdet*gauss->weight*(water_load);
			if(dt!=0.) scalar = scalar*dt;
			for(int i=0;i<numnodes;i++){
				pe->values[i]+=scalar*basis[i];
			}
		}
		else{
			/*if EPL is present and active input is there not here*/
			active_element_input->GetInputValue(&active_element);
			if(!active_element){
				water_input->GetInputValue(&water_load,gauss);
				scalar = Jdet*gauss->weight*(water_load);
				if(dt!=0.) scalar = scalar*dt;
				for(int i=0;i<numnodes;i++){
					pe->values[i]+=scalar*basis[i];
				}
			}
		}

		/*Transient and transfer terms*/
		if(dt!=0.){
			old_wh_input->GetInputValue(&water_head,gauss);
			sediment_storing = SedimentStoring(basalelement,gauss,sed_head_input,base_input);
			if(isefficientlayer){
				/*Dealing with the sediment part of the transfer term*/
				active_element_input->GetInputValue(&active_element);
				if(active_element){
					transfer=GetHydrologyPVectorTransfer(basalelement,gauss,epl_head_input);
				}
				else{
					transfer=0.0;
				}
				scalar = Jdet*gauss->weight*((water_head*sediment_storing)+(dt*transfer));
				for(int i=0;i<numnodes;i++)pe->values[i]+=scalar*basis[i];
			}
			else{
				scalar = Jdet*gauss->weight*(water_head*sediment_storing);
				for(int i=0;i<numnodes;i++)pe->values[i]+=scalar*basis[i];
			}
		}
	}
	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/

void HydrologyDCInefficientAnalysis::GetB(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
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

void HydrologyDCInefficientAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,SedimentHeadEnum);
}/*}}}*/

void HydrologyDCInefficientAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/

void HydrologyDCInefficientAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Intermediaries*/
	int						domaintype;
	Element*			basalelement=NULL;
	bool					converged;
	int*					doflist	=	NULL;

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


	/*Fetch number of nodes for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	basalelement->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values   = xNew<IssmDouble>(numnodes);
	IssmDouble* pressure = xNew<IssmDouble>(numnodes);
	IssmDouble* residual = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		values[i] =solution[doflist[i]];
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*If converged keep the residual in mind, also compute effective pressure*/
	basalelement->GetInputValue(&converged,ConvergedEnum);

	/*Get inputs*/
	if(converged){
		IssmDouble  penalty_factor,kmax,kappa,h_max;
		IssmDouble* thickness = xNew<IssmDouble>(numnodes);
		IssmDouble* base      = xNew<IssmDouble>(numnodes);

		basalelement->FindParam(&kmax,HydrologySedimentKmaxEnum);
		basalelement->FindParam(&penalty_factor,HydrologydcPenaltyFactorEnum);
		IssmDouble rho_freshwater = basalelement->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
		IssmDouble rho_ice        = basalelement->GetMaterialParameter(MaterialsRhoIceEnum);
		IssmDouble g              = basalelement->GetMaterialParameter(ConstantsGEnum);
		
		basalelement->GetInputListOnVertices(thickness,ThicknessEnum);
		basalelement->GetInputListOnVertices(base,BaseEnum);

		kappa=kmax*pow(10.,penalty_factor);

		for(int i=0;i<numnodes;i++){
			GetHydrologyDCInefficientHmax(&h_max,basalelement,basalelement->GetNode(i));
			if(values[i]>h_max) {
				residual[i] = kappa*(values[i]-h_max);
			}
			else{
				residual[i] = 0.;
			}
			//adding base in min to take into account heads under bed wich don't change N
			//pressure[i]=(rho_ice*g*thickness[i])-(rho_freshwater*g*(max((min(h_max,values[i])-base[i]),0.0)));
			pressure[i]=(rho_ice*g*thickness[i])-(rho_freshwater*g*(values[i]-base[i]));
		}
		xDelete<IssmDouble>(thickness);
		xDelete<IssmDouble>(base);
	}

	/*Add input to the element: */
	element->AddBasalInput(SedimentHeadEnum,values,P1Enum);
	element->AddBasalInput(SedimentHeadResidualEnum,residual,P1Enum);
	element->AddBasalInput(EffectivePressureEnum,pressure,P1Enum);

	/*Free ressources:*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(residual);
	xDelete<IssmDouble>(pressure);
	xDelete<int>(doflist);
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/

void HydrologyDCInefficientAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/

IssmDouble HydrologyDCInefficientAnalysis::SedimentStoring(Element* element,Gauss* gauss,Input* sed_head_input, Input* base_input){/*{{{*/
	IssmDouble sediment_storing;
	IssmDouble storing;
	IssmDouble base_elev,prestep_head,water_sheet;
	IssmDouble rho_freshwater           = element->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	IssmDouble g                        = element->GetMaterialParameter(ConstantsGEnum);
	IssmDouble sediment_porosity        = element->GetMaterialParameter(HydrologydcSedimentPorosityEnum);
	IssmDouble sediment_thickness       = element->GetMaterialParameter(HydrologydcSedimentThicknessEnum);
	IssmDouble sediment_compressibility = element->GetMaterialParameter(HydrologydcSedimentCompressibilityEnum);
	IssmDouble water_compressibility    = element->GetMaterialParameter(HydrologydcWaterCompressibilityEnum);
	base_input->GetInputValue(&base_elev,gauss);
	sed_head_input->GetInputValue(&prestep_head,gauss);
	water_sheet=max(0.0,(prestep_head-(base_elev-sediment_thickness)));

	storing=rho_freshwater*g*sediment_porosity*sediment_thickness*(water_compressibility+(sediment_compressibility/sediment_porosity));
	
	//Heavyside approximation (1/(1+exp(-2kx))) with k=25 centering at thickness minus 1%
	sediment_storing=(sediment_porosity*exp(-20.*(water_sheet-0.99*sediment_thickness))+storing)/(1+exp(-20.*(water_sheet-0.99*sediment_thickness)));

	return sediment_storing;
}/*}}}*/

IssmDouble HydrologyDCInefficientAnalysis::SedimentTransmitivity(Element* element,Gauss* gauss,Input* sed_head_input, Input* base_input,Input* SedTrans_input){/*{{{*/
	IssmDouble sediment_transmitivity;
	IssmDouble FullLayer_transmitivity;
	IssmDouble base_elev,prestep_head,water_sheet;
	IssmDouble sediment_thickness       = element->GetMaterialParameter(HydrologydcSedimentThicknessEnum);
	base_input->GetInputValue(&base_elev,gauss);
	sed_head_input->GetInputValue(&prestep_head,gauss);
	SedTrans_input->GetInputValue(&FullLayer_transmitivity,gauss);
	water_sheet=max(0.0,(prestep_head-(base_elev-sediment_thickness)));

	if (water_sheet<=sediment_thickness){
		sediment_transmitivity=FullLayer_transmitivity*water_sheet/sediment_thickness;
	}
	else{
		sediment_transmitivity=FullLayer_transmitivity;
	}
	return sediment_transmitivity;
}/*}}}*/

void  HydrologyDCInefficientAnalysis::GetHydrologyDCInefficientHmax(IssmDouble* ph_max,Element* element, Node* innode){/*{{{*/
	
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
		element->GetInputValue(&thickness,innode,ThicknessEnum);
		element->GetInputValue(&bed,innode,BaseEnum);
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

IssmDouble HydrologyDCInefficientAnalysis::GetHydrologyKMatrixTransfer(Element* element){/*{{{*/

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

IssmDouble HydrologyDCInefficientAnalysis::GetHydrologyPVectorTransfer(Element* element, Gauss* gauss, Input* epl_head_input){/*{{{*/

	int transfermethod;
	IssmDouble epl_head;
	IssmDouble leakage,transfer;
	element->FindParam(&transfermethod,HydrologydcTransferFlagEnum);

	/*Switch between the different transfer methods cases*/
	switch(transfermethod){
	case 0:
		/*Just keepping the transfer to zero*/
		transfer=0.0;
		break;
	case 1:
		_assert_(epl_head_input);
		epl_head_input->GetInputValue(&epl_head,gauss);
		element->FindParam(&leakage,HydrologydcLeakageFactorEnum);
		transfer=epl_head*leakage;
		break;
	default:
		_error_("no case higher than 1 for the Transfer method");
	}
	return transfer;
}/*}}}*/

void HydrologyDCInefficientAnalysis::ElementizeEplMask(FemModel* femmodel){/*{{{*/

	bool     element_active;
	Element* element=NULL;

	for(int i=0;i<femmodel->elements->Size();i++){
		element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
			Input* node_mask_input = element->GetInput(HydrologydcMaskEplactiveNodeEnum); _assert_(node_mask_input);
		
		if(node_mask_input->Max()>0.){
			element_active = true;
		}
		else{
			element_active = false;
		}
		element->AddInput(new BoolInput(HydrologydcMaskEplactiveEltEnum,element_active));
	}
}/*}}}*/
