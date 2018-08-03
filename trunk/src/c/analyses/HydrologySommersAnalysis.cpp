#include "./HydrologySommersAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Define 2 hardcoded parameters*/
#define OMEGA 0.001    // parameter controlling transition to nonlinear resistance in basal system (dimensionless)
#define NU    1.787e-6 //kinematic water viscosity m^2/s
#define CT    7.5e-8  // Clapeyron slope (K/Pa) 
#define CW    4.22e3   // specific heat capacity of water (J/kg/K)

/*Model processing*/
void HydrologySommersAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	if(hydrology_model!=HydrologysommersEnum) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spchead",HydrologySommersAnalysisEnum,P1Enum);

}/*}}}*/
void HydrologySommersAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Sommers?*/
	if(hydrology_model!=HydrologysommersEnum) return;

	/*Create discrete loads for Moulins*/
	CreateSingleNodeToElementConnectivity(iomodel);
	for(int i=0;i<iomodel->numberofvertices;i++){
		if (iomodel->domaintype!=Domain3DEnum){
			/*keep only this partition's nodes:*/
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(iomodel->loadcounter+i+1,i,iomodel,HydrologySommersAnalysisEnum));
			}
		}
		else if(reCast<int>(iomodel->Data("md.mesh.vertexonbase")[i])){
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(iomodel->loadcounter+i+1,i,iomodel,HydrologySommersAnalysisEnum));
			}	
		}
	}
	iomodel->DeleteData(1,"md.mesh.vertexonbase");

	/*Deal with Neumann BC*/
	int M,N;
	int *segments = NULL;
	iomodel->FetchData(&segments,&M,&N,"md.mesh.segments");

	/*Check that the size seem right*/
	_assert_(N==3); _assert_(M>=3);
	for(int i=0;i<M;i++){
		if(iomodel->my_elements[segments[i*3+2]-1]){
			loads->AddObject(new Neumannflux(iomodel->loadcounter+i+1,i,iomodel,segments,HydrologySommersAnalysisEnum));
		}
	}

	xDelete<int>(segments);

}/*}}}*/
void HydrologySommersAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Sommers?*/
	if(hydrology_model!=HydrologysommersEnum) return;

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,HydrologySommersAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  HydrologySommersAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void HydrologySommersAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
	int    hydrology_model,frictionlaw;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Sommers?*/
	if(hydrology_model!=HydrologysommersEnum) return;

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
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.mask.groundedice_levelset",MaskGroundediceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(elements,"md.basalforcings.geothermalflux",BasalforcingsGeothermalfluxEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.head",HydrologyHeadEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.gap_height",HydrologyGapHeightEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.englacial_input",HydrologyEnglacialInputEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.moulin_input",HydrologyMoulinInputEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.bump_spacing",HydrologyBumpSpacingEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.bump_height",HydrologyBumpHeightEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.reynolds",HydrologyReynoldsEnum);
	iomodel->FetchDataToInput(elements,"md.hydrology.neumannflux",HydrologyNeumannfluxEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vy",VyEnum);
	iomodel->FindConstant(&frictionlaw,"md.friction.law");

	/*Friction law variables*/
	switch(frictionlaw){
		case 1:
			iomodel->FetchDataToInput(elements,"md.friction.coefficient",FrictionCoefficientEnum);
			iomodel->FetchDataToInput(elements,"md.friction.p",FrictionPEnum);
			iomodel->FetchDataToInput(elements,"md.friction.q",FrictionQEnum);
			break;
		case 8:
			iomodel->FetchDataToInput(elements,"md.friction.coefficient",FrictionCoefficientEnum);
			break;
		default:
			_error_("Friction law "<< frictionlaw <<" not supported");
	}
}/*}}}*/
void HydrologySommersAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Sommers?*/
	if(hydrology_model!=HydrologysommersEnum) return;

	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));
	parameters->AddObject(iomodel->CopyConstantObject("md.friction.law",FrictionLawEnum));
   parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.relaxation",HydrologyRelaxationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.storage",HydrologyStorageEnum));
}/*}}}*/

/*Finite Element Analysis*/
void           HydrologySommersAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologySommersAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* HydrologySommersAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologySommersAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble Jdet;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/*Get conductivity from inputs*/
	IssmDouble conductivity = GetConductivity(element);

	/*Get englacial storage coefficient*/
	IssmDouble storage,dt;
	element->FindParam(&storage,HydrologyStorageEnum);
	element->FindParam(&dt,TimesteppingTimeStepEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(1);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += conductivity*gauss->weight*Jdet*(dbasis[0*numnodes+i]*dbasis[0*numnodes+j] + dbasis[1*numnodes+i]*dbasis[1*numnodes+j])
				  + gauss->weight*Jdet*storage/dt*basis[i]*basis[j];
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* HydrologySommersAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsFloating()) return NULL;

	/*Intermediaries */
	IssmDouble  Jdet,meltrate,G,dh[2],B,A,n;
	IssmDouble  gap,bed,thickness,head,ieb,head_old;
	IssmDouble  lr,br,vx,vy,beta;
	IssmDouble  alpha2,frictionheat;
   IssmDouble  PMPheat,dpressure_water[2],dbed[2];	
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	IssmDouble  latentheat      = element->GetMaterialParameter(MaterialsLatentheatEnum);
	IssmDouble  g               = element->GetMaterialParameter(ConstantsGEnum);
	IssmDouble  rho_ice         = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble  rho_water       = element->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	Input* geothermalflux_input = element->GetInput(BasalforcingsGeothermalfluxEnum);_assert_(geothermalflux_input);
	Input* head_input           = element->GetInput(HydrologyHeadEnum);              _assert_(head_input);
	Input* gap_input            = element->GetInput(HydrologyGapHeightEnum);         _assert_(gap_input);
	Input* thickness_input      = element->GetInput(ThicknessEnum);                  _assert_(thickness_input);
	Input* base_input           = element->GetInput(BaseEnum);                       _assert_(base_input);
	Input* B_input              = element->GetInput(MaterialsRheologyBEnum);         _assert_(B_input);
	Input* n_input              = element->GetInput(MaterialsRheologyNEnum);         _assert_(n_input);
	Input* englacial_input      = element->GetInput(HydrologyEnglacialInputEnum);    _assert_(englacial_input);
	Input* vx_input             = element->GetInput(VxEnum);                         _assert_(vx_input);
	Input* vy_input             = element->GetInput(VyEnum);                         _assert_(vy_input);
	Input* lr_input             = element->GetInput(HydrologyBumpSpacingEnum);       _assert_(lr_input);
	Input* br_input             = element->GetInput(HydrologyBumpHeightEnum);        _assert_(br_input);
   Input* headold_input        = element->GetInput(HydrologyHeadOldEnum);           _assert_(headold_input);

	/*Get conductivity from inputs*/
	IssmDouble conductivity = GetConductivity(element);

	/*Get englacial storage coefficient*/
	IssmDouble storage,dt;
   element->FindParam(&storage,HydrologyStorageEnum);
   element->FindParam(&dt,TimesteppingTimeStepEnum);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(element,2);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		geothermalflux_input->GetInputValue(&G,gauss);
		base_input->GetInputValue(&bed,gauss);
		base_input->GetInputDerivativeValue(&dbed[0],xyz_list,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		gap_input->GetInputValue(&gap,gauss);
		head_input->GetInputValue(&head,gauss);
		head_input->GetInputDerivativeValue(&dh[0],xyz_list,gauss);
		englacial_input->GetInputValue(&ieb,gauss);
		lr_input->GetInputValue(&lr,gauss);
		br_input->GetInputValue(&br,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
      headold_input->GetInputValue(&head_old,gauss);

		/*Get ice A parameter*/
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		A=pow(B,-n);
		
		/*Compute beta term*/
		if(gap<br)
		 beta = (br-gap)/lr;
		else
		 beta = 0.;

		/*Compute frictional heat flux*/
		friction->GetAlpha2(&alpha2,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		frictionheat=alpha2*(vx*vx+vy*vy);

		/*Get water and ice pressures*/
		IssmDouble pressure_ice   = rho_ice*g*thickness;    _assert_(pressure_ice>0.); 
		IssmDouble pressure_water = rho_water*g*(head-bed);
		if(pressure_water>pressure_ice) pressure_water = pressure_ice;

		/*Get water pressure from previous time step to use in lagged creep term*/
		IssmDouble pressure_water_old = rho_water*g*(head_old-bed);
		if(pressure_water_old>pressure_ice) pressure_water_old = pressure_ice;

		/*Compute change in sensible heat due to changes in pressure melting point*/
   	dpressure_water[0] = rho_water*g*(dh[0] - dbed[0]);
		dpressure_water[1] = rho_water*g*(dh[1] - dbed[1]);
		PMPheat=-CT*CW*conductivity*(dh[0]*dpressure_water[0]+dh[1]*dpressure_water[1]);

		meltrate = 1/latentheat*(G+frictionheat+rho_water*g*conductivity*(dh[0]*dh[0]+dh[1]*dh[1])-PMPheat);
		_assert_(meltrate>0.);
	
		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*
		 (
		  meltrate*(1/rho_water-1/rho_ice)
		  +A*pow(fabs(pressure_ice - pressure_water),n-1)*(pressure_ice - pressure_water)*gap
		  -beta*sqrt(vx*vx+vy*vy)
		  +ieb
		  +storage*head_old/dt
		  )*basis[i];     	
	}
	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete friction;
	delete gauss;
	return pe;
}/*}}}*/
void           HydrologySommersAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,HydrologyHeadEnum);
}/*}}}*/
void           HydrologySommersAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologySommersAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Intermediary*/
	IssmDouble dh[3];
	int* doflist = NULL;
	IssmDouble* xyz_list = NULL;

	/*Get gravity from parameters*/
	   IssmDouble  g = element->GetMaterialParameter(ConstantsGEnum);

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	element->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numnodes);

	/*Get thickness and base on nodes to apply cap on water head*/
   IssmDouble* eff_pressure = xNew<IssmDouble>(numnodes);
	IssmDouble* thickness = xNew<IssmDouble>(numnodes);
	IssmDouble* bed       = xNew<IssmDouble>(numnodes);
	IssmDouble  rho_ice   = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble  rho_water = element->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	element->GetInputListOnNodes(&thickness[0],ThicknessEnum);
	element->GetInputListOnNodes(&bed[0],BaseEnum);

	/*Get head from previous time-step and under-relaxation coefficient to use in under-relaxation for nonlinear convergence*/
   IssmDouble* head_old  = xNew<IssmDouble>(numnodes); 
	element->GetInputListOnNodes(&head_old[0],HydrologyHeadEnum);
   IssmDouble relaxation; 
	element->FindParam(&relaxation,HydrologyRelaxationEnum);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		values[i]=solution[doflist[i]];

		/*make sure that p_water<p_ice ->  h<rho_i H/rho_w + zb*/
		if(values[i]>rho_ice*thickness[i]/rho_water+bed[i]){
			values[i] = rho_ice*thickness[i]/rho_water+bed[i];
		}

		/*Make sure that negative pressure is not allowed*/
      if(values[i]<bed[i]){
			values[i] = bed[i];
		}

		/*Under-relaxation*/
	   values[i] = head_old[i] - relaxation*(head_old[i]-values[i]);

		/*Calculate effective pressure*/
		eff_pressure[i] = rho_ice*g*thickness[i] - rho_water*g*(values[i]-bed[i]);
	
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*Add input to the element: */
	element->AddInput(HydrologyHeadEnum,values,element->GetElementType());
   element->AddInput(EffectivePressureEnum,eff_pressure,P1Enum);

	/*Update reynolds number according to new solution*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* head_input = element->GetInput(HydrologyHeadEnum);_assert_(head_input);
	head_input->GetInputDerivativeAverageValue(&dh[0],xyz_list);
	IssmDouble conductivity = GetConductivity(element);

	IssmDouble reynolds = conductivity*sqrt(dh[0]*dh[0]+dh[1]*dh[1])/NU;
	element->AddInput(HydrologyReynoldsEnum,&reynolds,P0Enum);

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(thickness);
	xDelete<IssmDouble>(bed);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
	xDelete<IssmDouble>(eff_pressure);
   xDelete<IssmDouble>(head_old);
}/*}}}*/
void           HydrologySommersAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/

/*Additional methods*/
IssmDouble HydrologySommersAnalysis::GetConductivity(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble gap,reynolds;

	/*Get gravity from parameters*/
	IssmDouble  g = element->GetMaterialParameter(ConstantsGEnum);

	/*Get Reynolds and gap average values*/
	Input* reynolds_input = element->GetInput(HydrologyReynoldsEnum);  _assert_(reynolds_input);
	Input* gap_input      = element->GetInput(HydrologyGapHeightEnum); _assert_(gap_input);
	reynolds_input->GetInputAverage(&reynolds);
	gap_input->GetInputAverage(&gap);
	
	/*Compute conductivity*/
	IssmDouble conductivity = pow(gap,3)*g/(12.*NU*(1+OMEGA*reynolds));
	_assert_(conductivity>0);

	/*Clean up and return*/
	return conductivity;
}/*}}}*/
void HydrologySommersAnalysis::UpdateGapHeight(FemModel* femmodel){/*{{{*/


	for(int j=0;j<femmodel->elements->Size();j++){
		Element* element=(Element*)femmodel->elements->GetObjectByOffset(j);
		UpdateGapHeight(element);
	}

}/*}}}*/
void HydrologySommersAnalysis::UpdateGapHeight(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsFloating()) return;

	/*Intermediaries */
	IssmDouble newgap = 0.;
	IssmDouble  Jdet,meltrate,G,dh[2],B,A,n,dt;
	IssmDouble  gap,bed,thickness,head,ieb;
	IssmDouble  lr,br,vx,vy,beta;
	IssmDouble  alpha2,frictionheat;
	IssmDouble* xyz_list = NULL;
   IssmDouble  dpressure_water[2],dbed[2],PMPheat;
	IssmDouble q = 0.;

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	IssmDouble  latentheat      = element->GetMaterialParameter(MaterialsLatentheatEnum);
	IssmDouble  g               = element->GetMaterialParameter(ConstantsGEnum);
	IssmDouble  rho_ice         = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble  rho_water       = element->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	Input* geothermalflux_input = element->GetInput(BasalforcingsGeothermalfluxEnum);_assert_(geothermalflux_input);
	Input* head_input           = element->GetInput(HydrologyHeadEnum);              _assert_(head_input);
	Input* gap_input            = element->GetInput(HydrologyGapHeightEnum);         _assert_(gap_input);
	Input* thickness_input      = element->GetInput(ThicknessEnum);                  _assert_(thickness_input);
	Input* base_input           = element->GetInput(BaseEnum);                       _assert_(base_input);
	Input* B_input              = element->GetInput(MaterialsRheologyBEnum);         _assert_(B_input);
	Input* n_input              = element->GetInput(MaterialsRheologyNEnum);         _assert_(n_input);
	Input* englacial_input      = element->GetInput(HydrologyEnglacialInputEnum);    _assert_(englacial_input);
	Input* vx_input             = element->GetInput(VxEnum);                         _assert_(vx_input);
	Input* vy_input             = element->GetInput(VyEnum);                         _assert_(vy_input);
	Input* lr_input             = element->GetInput(HydrologyBumpSpacingEnum);       _assert_(lr_input);
	Input* br_input             = element->GetInput(HydrologyBumpHeightEnum);        _assert_(br_input);

	/*Get conductivity from inputs*/
	IssmDouble conductivity = GetConductivity(element);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(element,2);

	/*Keep track of weights*/
	IssmDouble totalweights=0.;

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);

		geothermalflux_input->GetInputValue(&G,gauss);
		base_input->GetInputValue(&bed,gauss);
		base_input->GetInputDerivativeValue(&dbed[0],xyz_list,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		gap_input->GetInputValue(&gap,gauss);
		head_input->GetInputValue(&head,gauss);
		head_input->GetInputDerivativeValue(&dh[0],xyz_list,gauss);
		englacial_input->GetInputValue(&ieb,gauss);
		lr_input->GetInputValue(&lr,gauss);
		br_input->GetInputValue(&br,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);

		/*Get ice A parameter*/
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		A=pow(B,-n);

		/*Compute beta term*/
		if(gap<br)
		 beta = (br-gap)/lr;
		else
		 beta = 0.;

		/*Compute frictional heat flux*/
		friction->GetAlpha2(&alpha2,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		frictionheat=alpha2*(vx*vx+vy*vy);

		/*Get water and ice pressures*/
		IssmDouble pressure_ice   = rho_ice*g*thickness;    _assert_(pressure_ice>0.); 
		IssmDouble pressure_water = rho_water*g*(head-bed);
		if(pressure_water>pressure_ice) pressure_water = pressure_ice;
      
      /* Compute change in sensible heat due to changes in pressure melting point*/
	   dpressure_water[0] = rho_water*g*(dh[0] - dbed[0]);
		dpressure_water[1] = rho_water*g*(dh[1] - dbed[1]);
		PMPheat=-CT*CW*conductivity*(dh[0]*dpressure_water[0]+dh[1]*dpressure_water[1]);
	
		meltrate = 1/latentheat*(G+frictionheat+rho_water*g*conductivity*(dh[0]*dh[0]+dh[1]*dh[1])-PMPheat);
		_assert_(meltrate>0.);

		newgap += gauss->weight*Jdet*(gap+dt*(
					meltrate/rho_ice
					-A*pow(fabs(pressure_ice-pressure_water),n-1)*(pressure_ice-pressure_water)*gap
					+beta*sqrt(vx*vx+vy*vy)
					));
		totalweights +=gauss->weight*Jdet;

		/* Compute basal water flux */
      q += gauss->weight*Jdet*(conductivity*sqrt(dh[0]*dh[0]+dh[1]*dh[1]));
	}

	/*Divide by connectivity*/
	newgap = newgap/totalweights;
	IssmDouble mingap = 1e-6;
	if(newgap<mingap) newgap=mingap;

	/*Limit gap height to grow to surface*/
	if(newgap>thickness)
	 newgap = thickness;
	 
	/*Add new gap as an input*/
	element->AddInput(HydrologyGapHeightEnum,&newgap,P0Enum);
 
	/*Divide by connectivity, add basal flux as an input*/
	q = q/totalweights;
	element->AddInput(HydrologyBasalFluxEnum,&q,P0Enum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	delete friction;
	delete gauss;
}/*}}}*/
