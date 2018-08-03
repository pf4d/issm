#include "./EnthalpyAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"
#include "../cores/cores.h"

/*Model processing*/
void EnthalpyAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Intermediary*/
	int        count;
	int        M,N;
	bool       spcpresent = false;
	int        finiteelement;
	IssmDouble heatcapacity;
	IssmDouble referencetemperature;

	/*Output*/
	IssmDouble *spcvector  = NULL;
	IssmDouble* times=NULL;
	IssmDouble* values=NULL;

	/*Fetch parameters: */
	iomodel->FindConstant(&heatcapacity,"md.materials.heatcapacity");
	iomodel->FindConstant(&referencetemperature,"md.constants.referencetemperature");
	iomodel->FindConstant(&finiteelement,"md.thermal.fe");

	/*return if 2d mesh*/
	if(iomodel->domaintype==Domain2DhorizontalEnum) return;

	/*Fetch data: */
	iomodel->FetchData(&spcvector,&M,&N,"md.thermal.spctemperature");

	/*Convert spcs from temperatures to enthalpy*/
	_assert_(N>0); _assert_(M>=iomodel->numberofvertices);
	for(int i=0;i<iomodel->numberofvertices;i++){
		for(int j=0;j<N;j++){
			spcvector[i*N+j] = heatcapacity*(spcvector[i*N+j]-referencetemperature);
		}
	}

	/*Specific case for PDD, we want the constaints to be updated by the PDD scheme itself*/
	bool isdynamic = false;
	if (iomodel->solution_enum==TransientSolutionEnum){
		int smb_model;
		iomodel->FindConstant(&smb_model,"md.smb.model");
		if(smb_model==SMBpddEnum)     isdynamic=true;
		if(smb_model==SMBd18opddEnum) isdynamic=true;
	}

	if(isdynamic){
		IoModelToDynamicConstraintsx(constraints,iomodel,spcvector,M,N,EnthalpyAnalysisEnum,finiteelement);
	}
	else{
		IoModelToConstraintsx(constraints,iomodel,spcvector,M,N,EnthalpyAnalysisEnum,finiteelement);
	}

	/*Free ressources:*/
	iomodel->DeleteData(spcvector,"md.thermal.spctemperature");
	xDelete<IssmDouble>(times);
	xDelete<IssmDouble>(values);
}/*}}}*/
void EnthalpyAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*No loads */
}/*}}}*/
void EnthalpyAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.thermal.fe");

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,EnthalpyAnalysisEnum,finiteelement);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  EnthalpyAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void EnthalpyAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	bool dakota_analysis,ismovingfront,isenthalpy;
	int frictionlaw,basalforcing_model,materialstype;
	int FrictionCoupling;
	
	/*Now, is the model 3d? otherwise, do nothing: */
	if(iomodel->domaintype==Domain2DhorizontalEnum)return;

	/*Is enthalpy requested?*/
	iomodel->FindConstant(&isenthalpy,"md.thermal.isenthalpy");
	if(!isenthalpy) return;

	/*Fetch data needed: */
	iomodel->FetchData(3,"md.initialization.temperature","md.initialization.waterfraction","md.initialization.pressure");

	/*Finite element type*/
	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.thermal.fe");

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");
	iomodel->FindConstant(&ismovingfront,"md.transient.ismovingfront");
	iomodel->FindConstant(&frictionlaw,"md.friction.law");
	iomodel->FindConstant(&materialstype,"md.materials.type");

	iomodel->FetchDataToInput(elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(elements,"md.slr.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.mask.groundedice_levelset",MaskGroundediceLevelsetEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(elements,"md.initialization.pressure",PressureEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.temperature",TemperatureEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.waterfraction",WaterfractionEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.enthalpy",EnthalpyEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.watercolumn",WatercolumnEnum);
	iomodel->FetchDataToInput(elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vy",VyEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vz",VzEnum);
	InputUpdateFromConstantx(elements,0.,VxMeshEnum);
	InputUpdateFromConstantx(elements,0.,VyMeshEnum);
	InputUpdateFromConstantx(elements,0.,VzMeshEnum);
	if(ismovingfront){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum); // required for updating active nodes
	}

	/*Basal forcings variables*/
	iomodel->FindConstant(&basalforcing_model,"md.basalforcings.model");
	switch(basalforcing_model){
		case MantlePlumeGeothermalFluxEnum:
			break;
		default:
			iomodel->FetchDataToInput(elements,"md.basalforcings.geothermalflux",BasalforcingsGeothermalfluxEnum);
			break;
	}

	/*Rheology type*/
	iomodel->FetchDataToInput(elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
	switch(materialstype){
		case MatenhancediceEnum:
			iomodel->FetchDataToInput(elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			iomodel->FetchDataToInput(elements,"md.materials.rheology_E",MaterialsRheologyEEnum);
			break;
		case MatdamageiceEnum:
			iomodel->FetchDataToInput(elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			break;
		case MatestarEnum:
			iomodel->FetchDataToInput(elements,"md.materials.rheology_Ec",MaterialsRheologyEcEnum);
			iomodel->FetchDataToInput(elements,"md.materials.rheology_Es",MaterialsRheologyEsEnum);
			break;
		case MaticeEnum:
			iomodel->FetchDataToInput(elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			break;
		default:
			_error_("not supported");
	}
	
	/*Friction law variables*/
	switch(frictionlaw){
		case 1:
			iomodel->FetchDataToInput(elements,"md.friction.coefficient",FrictionCoefficientEnum);
			iomodel->FetchDataToInput(elements,"md.friction.p",FrictionPEnum);
			iomodel->FetchDataToInput(elements,"md.friction.q",FrictionQEnum);
			break;
		case 2:
			iomodel->FetchDataToInput(elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(elements,"md.friction.m",FrictionMEnum);
			break;
		case 3:
			iomodel->FindConstant(&FrictionCoupling,"md.friction.coupling");
			iomodel->FetchDataToInput(elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(elements,"md.friction.As",FrictionAsEnum);
			iomodel->FetchDataToInput(elements,"md.friction.q",FrictionQEnum);
			if (FrictionCoupling==0){
				iomodel->FetchDataToInput(elements,"md.friction.effective_pressure",FrictionEffectivePressureEnum);
			}
			break;
		case 4:
			iomodel->FetchDataToInput(elements,"md.friction.coefficient",FrictionCoefficientEnum);
			iomodel->FetchDataToInput(elements,"md.friction.p",FrictionPEnum);
			iomodel->FetchDataToInput(elements,"md.friction.q",FrictionQEnum);
			iomodel->FetchDataToInput(elements,"md.initialization.pressure",PressureEnum);
			iomodel->FetchDataToInput(elements,"md.initialization.temperature",TemperatureEnum);
			break;
		case 5:
			iomodel->FetchDataToInput(elements,"md.friction.coefficient",FrictionCoefficientEnum);
			iomodel->FetchDataToInput(elements,"md.friction.p",FrictionPEnum);
			iomodel->FetchDataToInput(elements,"md.friction.q",FrictionQEnum);
			iomodel->FetchDataToInput(elements,"md.friction.water_layer",FrictionWaterLayerEnum);
			break;
		case 6:
			iomodel->FetchDataToInput(elements,"md.friction.C",FrictionCEnum);
			iomodel->FetchDataToInput(elements,"md.friction.m",FrictionMEnum);
			iomodel->FetchDataToInput(elements,"md.initialization.pressure",PressureEnum);
			iomodel->FetchDataToInput(elements,"md.initialization.temperature",TemperatureEnum);
			break;
		case 9:
			iomodel->FetchDataToInput(elements,"md.friction.coefficient",FrictionCoefficientEnum);
			iomodel->FetchDataToInput(elements,"md.friction.pressure_adjusted_temperature",FrictionPressureAdjustedTemperatureEnum);
			InputUpdateFromConstantx(elements,1.,FrictionPEnum);
			InputUpdateFromConstantx(elements,1.,FrictionQEnum);
			break;
		default:
			_error_("not supported");
	}

	/*Free data: */
	iomodel->DeleteData(3,"md.initialization.temperature","md.initialization.waterfraction","md.initialization.pressure");

}/*}}}*/
void EnthalpyAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;

	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.stabilization",ThermalStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.maxiter",ThermalMaxiterEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.reltol",ThermalReltolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.isenthalpy",ThermalIsenthalpyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.thermal.isdynamicbasalspc",ThermalIsdynamicbasalspcEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.friction.law",FrictionLawEnum));

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.thermal.requested_outputs");
	parameters->AddObject(new IntParam(ThermalNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(ThermalRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.thermal.requested_outputs");

	/*Deal with friction parameters*/
	int frictionlaw;
	iomodel->FindConstant(&frictionlaw,"md.friction.law");
	if(frictionlaw==4 || frictionlaw==6) parameters->AddObject(iomodel->CopyConstantObject("md.friction.gamma",FrictionGammaEnum));
	if(frictionlaw==3) parameters->AddObject(iomodel->CopyConstantObject("md.friction.coupling",FrictionCouplingEnum));
	if(frictionlaw==9) parameters->AddObject(iomodel->CopyConstantObject("md.friction.gamma",FrictionGammaEnum));
}/*}}}*/

/*Finite Element Analysis*/
void           EnthalpyAnalysis::ApplyBasalConstraints(IssmDouble* serial_spc,Element* element){/*{{{*/

	/* Do not check if ice in element, this may lead to inconsistencies between cpu partitions */
	/* Only update constraints at the base. */
	if(!(element->IsOnBase())) return;

	/*Intermediary*/
	bool        isdynamicbasalspc;
	int         numindices;
	int        *indices = NULL;
	Node*       node = NULL;
	IssmDouble	pressure;

	/*Check wether dynamic basal boundary conditions are activated */
	element->FindParam(&isdynamicbasalspc,ThermalIsdynamicbasalspcEnum);
	if(!isdynamicbasalspc) return;

	/*Get parameters and inputs: */
	Input* pressure_input		 = element->GetInput(PressureEnum);							 _assert_(pressure_input);

	/*Fetch indices of basal & surface nodes for this finite element*/
	Penta *penta =  (Penta *) element; // TODO: add Basal-/SurfaceNodeIndices to element.h, and change this to Element*
	penta->BasalNodeIndices(&numindices,&indices,element->GetElementType());

	GaussPenta* gauss=new GaussPenta();
	for(int i=0;i<numindices;i++){
		gauss->GaussNode(element->GetElementType(),indices[i]);

		pressure_input->GetInputValue(&pressure,gauss);

		/*apply or release spc*/
		node=element->GetNode(indices[i]);
		if(!node->IsActive()) continue;
		if(serial_spc[node->Sid()]==1.){
			pressure_input->GetInputValue(&pressure, gauss);
			node->ApplyConstraint(0,PureIceEnthalpy(element,pressure));
		}
		else {
			node->DofInFSet(0);
		}
	}

	/*Free ressources:*/
	xDelete<int>(indices);
	delete gauss;
}/*}}}*/
void           EnthalpyAnalysis::ComputeBasalMeltingrate(FemModel* femmodel){/*{{{*/
	/*Compute basal melting rates: */
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ComputeBasalMeltingrate(element);
	}
}/*}}}*/
void           EnthalpyAnalysis::ComputeBasalMeltingrate(Element* element){/*{{{*/
	/*Calculate the basal melt rates of the enthalpy model after Aschwanden 2012*/
	/* melting rate is positive when melting, negative when refreezing*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return;

	/* Only compute melt rates at the base of grounded ice*/
	if(!element->IsOnBase() || element->IsFloating()) return;

	/* Intermediaries */
	bool			converged;
	const int   dim=3;
	int         i,is,state;
	int			nodedown,nodeup,numnodes,numsegments;
	int			enthalpy_enum;
	IssmDouble  vec_heatflux[dim],normal_base[dim],d1enthalpy[dim],d1pressure[dim];
	IssmDouble  basalfriction,alpha2,geothermalflux,heatflux;
	IssmDouble  dt,yts;
	IssmDouble  melting_overshoot,lambda;
	IssmDouble  vx,vy,vz;
	IssmDouble *xyz_list      = NULL;
	IssmDouble *xyz_list_base = NULL;
	int        *pairindices   = NULL;

	/*Fetch parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->GetInputValue(&converged,ConvergedEnum);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&yts, ConstantsYtsEnum);

	if(dt==0. && !converged) enthalpy_enum=EnthalpyPicardEnum;
	else enthalpy_enum=EnthalpyEnum;

	IssmDouble latentheat = element->GetMaterialParameter(MaterialsLatentheatEnum);
	IssmDouble rho_ice    = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble rho_water  = element->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	IssmDouble beta		 = element->GetMaterialParameter(MaterialsBetaEnum);
	IssmDouble kappa		 = EnthalpyDiffusionParameterVolume(element,enthalpy_enum);     _assert_(kappa>=0.);
	IssmDouble kappa_mix;

	/*retrieve inputs*/
	Input* enthalpy_input         = element->GetInput(enthalpy_enum);                    _assert_(enthalpy_input);
	Input* pressure_input			= element->GetInput(PressureEnum);							 _assert_(pressure_input);
	Input* geothermalflux_input   = element->GetInput(BasalforcingsGeothermalfluxEnum); _assert_(geothermalflux_input);
	Input* vx_input               = element->GetInput(VxEnum);                          _assert_(vx_input);
	Input* vy_input               = element->GetInput(VyEnum);                          _assert_(vy_input);
	Input* vz_input               = element->GetInput(VzEnum);                          _assert_(vz_input);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(element,dim);

	/******** MELTING RATES  ************************************//*{{{*/
	element->NormalBase(&normal_base[0],xyz_list_base);
	element->VerticalSegmentIndicesBase(&pairindices,&numsegments);
	IssmDouble* meltingrate_enthalpy = xNew<IssmDouble>(numsegments);
	IssmDouble* heating = xNew<IssmDouble>(numsegments);	

	numnodes=element->GetNumberOfNodes();
	IssmDouble* enthalpies = xNew<IssmDouble>(numnodes);
	IssmDouble* pressures = xNew<IssmDouble>(numnodes);
	IssmDouble* watercolumns = xNew<IssmDouble>(numnodes);
	IssmDouble* basalmeltingrates = xNew<IssmDouble>(numnodes);
	element->GetInputListOnNodes(enthalpies,enthalpy_enum);
	element->GetInputListOnNodes(pressures,PressureEnum);
	element->GetInputListOnNodes(watercolumns,WatercolumnEnum);
	element->GetInputListOnNodes(basalmeltingrates,BasalforcingsGroundediceMeltingRateEnum);

	Gauss* gauss=element->NewGauss();
	for(is=0;is<numsegments;is++){
		nodedown = pairindices[is*2+0];
		nodeup   = pairindices[is*2+1];
		gauss->GaussNode(element->GetElementType(),nodedown);

		state=GetThermalBasalCondition(element, enthalpies[nodedown], enthalpies[nodeup], pressures[nodedown], pressures[nodeup], watercolumns[nodedown], basalmeltingrates[nodedown]);
		switch (state) {
			case 0:
				// cold, dry base: apply basal surface forcing
				for(i=0;i<3;i++) vec_heatflux[i]=0.;
				break;
			case 1: case 2: case 3: 
				// case 1 : cold, wet base: keep at pressure melting point 
				// case 2: temperate, thin refreezing base: release spc
				// case 3: temperate, thin melting base: set spc
				enthalpy_input->GetInputDerivativeValue(&d1enthalpy[0],xyz_list,gauss);
				for(i=0;i<3;i++) vec_heatflux[i]=-kappa*d1enthalpy[i];
				break;
			case 4:
				// temperate, thick melting base: set grad H*n=0
				kappa_mix=GetWetIceConductivity(element, enthalpies[nodedown], pressures[nodedown]);
				pressure_input->GetInputDerivativeValue(&d1pressure[0],xyz_list,gauss);
				for(i=0;i<3;i++) vec_heatflux[i]=kappa_mix*beta*d1pressure[i];
				break;
			default:
				_printf0_("	unknown thermal basal state found!");
		}
		if(state==0) meltingrate_enthalpy[is]=0.;
		else{
			/*heat flux along normal*/
			heatflux=0.;
			for(i=0;i<3;i++) heatflux+=(vec_heatflux[i])*normal_base[i];

			/*basal friction*/
			friction->GetAlpha2(&alpha2,gauss);
			vx_input->GetInputValue(&vx,gauss);		vy_input->GetInputValue(&vy,gauss);		vz_input->GetInputValue(&vz,gauss);
			basalfriction=alpha2*(vx*vx + vy*vy + vz*vz);
			geothermalflux_input->GetInputValue(&geothermalflux,gauss);
			/* -Mb= Fb-(q-q_geo)/((1-w)*L*rho), and (1-w)*rho=rho_ice, cf Aschwanden 2012, eqs.1, 2, 66*/
			heating[is]=(heatflux+basalfriction+geothermalflux);
			meltingrate_enthalpy[is]=heating[is]/(latentheat*rho_ice); // m/s water equivalent
		}
	}/*}}}*/

	/******** UPDATE MELTINGRATES AND WATERCOLUMN **************//*{{{*/
	for(is=0;is<numsegments;is++){
		nodedown = pairindices[is*2+0];
		nodeup   = pairindices[is*2+1];
		if(dt!=0.){
			if(watercolumns[nodedown]+meltingrate_enthalpy[is]*dt<0.){	// prevent too much freeze on			
				lambda = -watercolumns[nodedown]/(dt*meltingrate_enthalpy[is]); _assert_(lambda>=0.); _assert_(lambda<1.);
				watercolumns[nodedown]=0.;
				basalmeltingrates[nodedown]=lambda*meltingrate_enthalpy[is]; // restrict freeze on only to size of watercolumn
				enthalpies[nodedown]+=(1.-lambda)*dt/yts*meltingrate_enthalpy[is]*latentheat*rho_ice; // use rest of energy to cool down base: dE=L*m, m=(1-lambda)*meltingrate*rho_ice
			}
			else{
				basalmeltingrates[nodedown]=meltingrate_enthalpy[is];
				watercolumns[nodedown]+=dt*meltingrate_enthalpy[is]; 
			}
		}
		else{
			basalmeltingrates[nodedown]=meltingrate_enthalpy[is];
			if(watercolumns[nodedown]+meltingrate_enthalpy[is]<0.)
				watercolumns[nodedown]=0.;
			else
				watercolumns[nodedown]+=meltingrate_enthalpy[is];
		}	
		basalmeltingrates[nodedown]*=rho_water/rho_ice; // convert meltingrate from water to ice equivalent
		_assert_(watercolumns[nodedown]>=0.);
	}/*}}}*/

	/*feed updated variables back into model*/
	if(dt!=0.){
		element->AddInput(enthalpy_enum,enthalpies,element->GetElementType()); 
		element->AddInput(WatercolumnEnum,watercolumns,element->GetElementType());
	}
	element->AddInput(BasalforcingsGroundediceMeltingRateEnum,basalmeltingrates,element->GetElementType());

	/*Clean up and return*/
	delete gauss;
	delete friction;
	xDelete<int>(pairindices);
	xDelete<IssmDouble>(enthalpies);
	xDelete<IssmDouble>(pressures);
	xDelete<IssmDouble>(watercolumns);
	xDelete<IssmDouble>(basalmeltingrates);
	xDelete<IssmDouble>(meltingrate_enthalpy);
	xDelete<IssmDouble>(heating);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_base);
}/*}}}*/
void           EnthalpyAnalysis::Core(FemModel* femmodel){/*{{{*/

	IssmDouble dt;
	bool isdynamicbasalspc;

	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	femmodel->parameters->FindParam(&isdynamicbasalspc,ThermalIsdynamicbasalspcEnum);

	if(VerboseSolution()) _printf0_("   computing enthalpy\n");
	femmodel->SetCurrentConfiguration(EnthalpyAnalysisEnum);
	if((dt>0.) && isdynamicbasalspc)	UpdateBasalConstraints(femmodel);
	solutionsequence_thermal_nonlinear(femmodel);

	/*transfer enthalpy to enthalpy picard for the next step: */
	InputDuplicatex(femmodel,EnthalpyEnum,EnthalpyPicardEnum);

	PostProcessing(femmodel);

}/*}}}*/
ElementVector* EnthalpyAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* EnthalpyAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* EnthalpyAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixVolume(element);
	ElementMatrix* Ke2=CreateKMatrixShelf(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	return Ke;
}/*}}}*/
ElementMatrix* EnthalpyAnalysis::CreateKMatrixVolume(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         stabilization;
	IssmDouble  Jdet,dt,u,v,w,um,vm,wm,vel;
	IssmDouble  h,hx,hy,hz,vx,vy,vz;
	IssmDouble  tau_parameter,diameter;
	IssmDouble  D_scalar;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke       = element->NewElementMatrix();
	IssmDouble*    basis    = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis   = xNew<IssmDouble>(3*numnodes);
	IssmDouble*    B        = xNew<IssmDouble>(3*numnodes);
	IssmDouble*    Bprime   = xNew<IssmDouble>(3*numnodes);
	IssmDouble     D[3][3]  = {0.};
	IssmDouble     K[3][3];

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&stabilization,ThermalStabilizationEnum);
	IssmDouble  rho_water           = element->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	IssmDouble  rho_ice             = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble  gravity             = element->GetMaterialParameter(ConstantsGEnum);
	IssmDouble  heatcapacity        = element->GetMaterialParameter(MaterialsHeatcapacityEnum);
	IssmDouble  thermalconductivity = element->GetMaterialParameter(MaterialsThermalconductivityEnum);
	Input* vx_input  = element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input  = element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input  = element->GetInput(VzEnum);     _assert_(vz_input);
	Input* vxm_input = element->GetInput(VxMeshEnum); _assert_(vxm_input);
	Input* vym_input = element->GetInput(VyMeshEnum); _assert_(vym_input);
	Input* vzm_input = element->GetInput(VzMeshEnum); _assert_(vzm_input);
	if(stabilization==2) diameter=element->MinEdgeLength(xyz_list);

	/*Enthalpy diffusion parameter*/
	IssmDouble kappa=this->EnthalpyDiffusionParameterVolume(element,EnthalpyPicardEnum); _assert_(kappa>=0.);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		D_scalar=gauss->weight*Jdet;
		if(dt!=0.) D_scalar=D_scalar*dt;

		/*Conduction: */
		GetBConduct(B,element,xyz_list,gauss); 
		D[0][0]=D_scalar*kappa/rho_ice;
		D[1][1]=D_scalar*kappa/rho_ice;
		D[2][2]=D_scalar*kappa/rho_ice;
		TripleMultiply(B,3,numnodes,1,
					&D[0][0],3,3,0,
					B,3,numnodes,0,
					&Ke->values[0],1);

		/*Advection: */
		GetBAdvec(B,element,xyz_list,gauss); 
		GetBAdvecprime(Bprime,element,xyz_list,gauss); 
		vx_input->GetInputValue(&u,gauss); vxm_input->GetInputValue(&um,gauss); vx=u-um;
		vy_input->GetInputValue(&v,gauss); vym_input->GetInputValue(&vm,gauss); vy=v-vm;
		vz_input->GetInputValue(&w,gauss); vzm_input->GetInputValue(&wm,gauss); vz=w-wm;
		D[0][0]=D_scalar*vx;
		D[1][1]=D_scalar*vy;
		D[2][2]=D_scalar*vz;
		TripleMultiply(B,3,numnodes,1,
					&D[0][0],3,3,0,
					Bprime,3,numnodes,0,
					&Ke->values[0],1);

		/*Transient: */
		if(dt!=0.){
			D_scalar=gauss->weight*Jdet;
			element->NodalFunctions(basis,gauss);
			TripleMultiply(basis,numnodes,1,0,
						&D_scalar,1,1,0,
						basis,1,numnodes,0,
						&Ke->values[0],1);
			D_scalar=D_scalar*dt;
		}

		/*Artificial diffusivity*/
		if(stabilization==1){
			element->ElementSizes(&hx,&hy,&hz);
			vel=sqrt(vx*vx + vy*vy + vz*vz)+1.e-14;
			h=sqrt( pow(hx*vx/vel,2) + pow(hy*vy/vel,2) + pow(hz*vz/vel,2));
			K[0][0]=h/(2.*vel)*fabs(vx*vx);  K[0][1]=h/(2.*vel)*fabs(vx*vy); K[0][2]=h/(2.*vel)*fabs(vx*vz);
			K[1][0]=h/(2.*vel)*fabs(vy*vx);  K[1][1]=h/(2.*vel)*fabs(vy*vy); K[1][2]=h/(2.*vel)*fabs(vy*vz);
			K[2][0]=h/(2.*vel)*fabs(vz*vx);  K[2][1]=h/(2.*vel)*fabs(vz*vy); K[2][2]=h/(2.*vel)*fabs(vz*vz);
			for(int i=0;i<3;i++) for(int j=0;j<3;j++) K[i][j] = D_scalar*K[i][j];

			GetBAdvecprime(Bprime,element,xyz_list,gauss); 
			TripleMultiply(Bprime,3,numnodes,1,
						&K[0][0],3,3,0,
						Bprime,3,numnodes,0,
						&Ke->values[0],1);
		}
		else if(stabilization==2){
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
			tau_parameter=element->StabilizationParameter(u-um,v-vm,w-wm,diameter,kappa/rho_ice);
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j]+=tau_parameter*D_scalar*
					  ((u-um)*dbasis[0*numnodes+i]+(v-vm)*dbasis[1*numnodes+i]+(w-wm)*dbasis[2*numnodes+i])*((u-um)*dbasis[0*numnodes+j]+(v-vm)*dbasis[1*numnodes+j]+(w-wm)*dbasis[2*numnodes+j]);
				}
			}
			if(dt!=0.){
				D_scalar=gauss->weight*Jdet;
				for(int i=0;i<numnodes;i++){
					for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j]+=tau_parameter*D_scalar*basis[j]*((u-um)*dbasis[0*numnodes+i]+(v-vm)*dbasis[1*numnodes+i]+(w-wm)*dbasis[2*numnodes+i]);
					}
				}
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(Bprime);
	delete gauss;
	return Ke;
}/*}}}*/
ElementMatrix* EnthalpyAnalysis::CreateKMatrixShelf(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Initialize Element matrix and return if necessary*/
	if(!element->IsOnBase() || !element->IsFloating()) return NULL;

	/*Intermediaries*/
	IssmDouble  dt,Jdet,D;
	IssmDouble *xyz_list_base = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize vectors*/
	ElementMatrix* Ke    = element->NewElementMatrix();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	IssmDouble  gravity             = element->GetMaterialParameter(ConstantsGEnum);
	IssmDouble  rho_water           = element->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	IssmDouble  rho_ice             = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble  heatcapacity        = element->GetMaterialParameter(MaterialsHeatcapacityEnum);
	IssmDouble  mixed_layer_capacity= element->GetMaterialParameter(MaterialsMixedLayerCapacityEnum);
	IssmDouble  thermal_exchange_vel= element->GetMaterialParameter(MaterialsThermalExchangeVelocityEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(4);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctions(basis,gauss);

		D=gauss->weight*Jdet*rho_water*mixed_layer_capacity*thermal_exchange_vel/(heatcapacity*rho_ice);
		if(reCast<bool,IssmDouble>(dt)) D=dt*D;
		TripleMultiply(basis,numnodes,1,0,
					&D,1,1,0,
					basis,1,numnodes,0,
					&Ke->values[0],1);

	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list_base);
	return Ke;
}/*}}}*/
ElementVector* EnthalpyAnalysis::CreatePVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*compute all load vectors for this element*/
	ElementVector* pe1=CreatePVectorVolume(element);
	ElementVector* pe2=CreatePVectorSheet(element);
	ElementVector* pe3=CreatePVectorShelf(element);
	ElementVector* pe =new ElementVector(pe1,pe2,pe3);

	/*clean-up and return*/
	delete pe1;
	delete pe2;
	delete pe3;
	return pe;
}/*}}}*/
ElementVector* EnthalpyAnalysis::CreatePVectorVolume(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int         i, stabilization;
	IssmDouble  Jdet,phi,dt;
	IssmDouble  enthalpy, Hpmp;
	IssmDouble  enthalpypicard, d1enthalpypicard[3];
	IssmDouble  pressure, d1pressure[3], d2pressure;
	IssmDouble  waterfractionpicard;
	IssmDouble  kappa,tau_parameter,diameter,kappa_w;
	IssmDouble  u,v,w;
	IssmDouble  scalar_def, scalar_sens ,scalar_transient;
	IssmDouble* xyz_list = NULL;
	IssmDouble  d1H_d1P, d1P2;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = element->GetNumberOfNodes();
	int numvertices = element->GetNumberOfVertices();

	/*Initialize Element vector*/
	ElementVector* pe             = element->NewElementVector();
	IssmDouble*    basis          = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis         = xNew<IssmDouble>(3*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	IssmDouble  rho_ice             = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble  heatcapacity        = element->GetMaterialParameter(MaterialsHeatcapacityEnum);
	IssmDouble  thermalconductivity = element->GetMaterialParameter(MaterialsThermalconductivityEnum);
	IssmDouble  temperateiceconductivity = element->GetMaterialParameter(MaterialsTemperateiceconductivityEnum);
	IssmDouble  beta                = element->GetMaterialParameter(MaterialsBetaEnum);
	IssmDouble  latentheat          = element->GetMaterialParameter(MaterialsLatentheatEnum);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&stabilization,ThermalStabilizationEnum);
	Input* vx_input=element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input=element->GetInput(VzEnum); _assert_(vz_input);
	Input* enthalpypicard_input=element->GetInput(EnthalpyPicardEnum); _assert_(enthalpypicard_input);
	Input* pressure_input=element->GetInput(PressureEnum); _assert_(pressure_input);
	Input* enthalpy_input=NULL;
	if(reCast<bool,IssmDouble>(dt)){enthalpy_input = element->GetInput(EnthalpyEnum); _assert_(enthalpy_input);}
	if(stabilization==2){
		diameter=element->MinEdgeLength(xyz_list);
		kappa=this->EnthalpyDiffusionParameterVolume(element,EnthalpyPicardEnum); _assert_(kappa>=0.);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(4);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		
		/*viscous dissipation*/
		element->ViscousHeating(&phi,xyz_list,gauss,vx_input,vy_input,vz_input);

		scalar_def=phi/rho_ice*Jdet*gauss->weight;
		if(dt!=0.) scalar_def=scalar_def*dt;

		for(i=0;i<numnodes;i++) pe->values[i]+=scalar_def*basis[i];

		/*sensible heat flux in temperate ice*/
		enthalpypicard_input->GetInputValue(&enthalpypicard,gauss);
		pressure_input->GetInputValue(&pressure,gauss);
		Hpmp=this->PureIceEnthalpy(element, pressure);

		if(enthalpypicard>=Hpmp){
			enthalpypicard_input->GetInputDerivativeValue(&d1enthalpypicard[0],xyz_list,gauss);
			pressure_input->GetInputDerivativeValue(&d1pressure[0],xyz_list,gauss);
			d2pressure=0.; // for linear elements, 2nd derivative is zero
			
			d1H_d1P=0.;
			for(i=0;i<3;i++) d1H_d1P+=d1enthalpypicard[i]*d1pressure[i];
			d1P2=0.;
			for(i=0;i<3;i++) d1P2+=pow(d1pressure[i],2.);

			scalar_sens=-beta*((temperateiceconductivity - thermalconductivity)/latentheat*(d1H_d1P + beta*heatcapacity*d1P2))/rho_ice;
			if(dt!=0.) scalar_sens=scalar_sens*dt;
			for(i=0;i<numnodes;i++) pe->values[i]+=scalar_sens*basis[i];
		}		

		/* Build transient now */
		if(reCast<bool,IssmDouble>(dt)){
			enthalpy_input->GetInputValue(&enthalpy, gauss);
			scalar_transient=enthalpy*Jdet*gauss->weight;
			for(i=0;i<numnodes;i++) pe->values[i]+=scalar_transient*basis[i];
		}

		if(stabilization==2){
			element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

			vx_input->GetInputValue(&u,gauss);
			vy_input->GetInputValue(&v,gauss);
			vz_input->GetInputValue(&w,gauss);
			tau_parameter=element->StabilizationParameter(u,v,w,diameter,kappa/rho_ice);

			for(i=0;i<numnodes;i++) pe->values[i]+=tau_parameter*scalar_def*(u*dbasis[0*numnodes+i]+v*dbasis[1*numnodes+i]+w*dbasis[2*numnodes+i]);

			if(dt!=0.){
				for(i=0;i<numnodes;i++) pe->values[i]+=tau_parameter*scalar_transient*(u*dbasis[0*numnodes+i]+v*dbasis[1*numnodes+i]+w*dbasis[2*numnodes+i]);
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return pe;

}/*}}}*/
ElementVector* EnthalpyAnalysis::CreatePVectorSheet(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/* implementation of the basal condition decision chart of Aschwanden 2012, Fig.5 */
	if(!element->IsOnBase() || element->IsFloating()) return NULL;

	bool converged, isdynamicbasalspc;
	int i, state;
	int enthalpy_enum;
	IssmDouble  dt,Jdet,scalar;
	IssmDouble	enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate;
	IssmDouble	vx,vy,vz;
	IssmDouble  alpha2,basalfriction,geothermalflux,heatflux;
	IssmDouble *xyz_list_base = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&isdynamicbasalspc,ThermalIsdynamicbasalspcEnum);
	element->GetInputValue(&converged,ConvergedEnum);
	if(dt==0. && !converged) enthalpy_enum=EnthalpyPicardEnum; // use enthalpy from last iteration
	else enthalpy_enum=EnthalpyEnum; // use enthalpy from last time step
	Input* vx_input             = element->GetInput(VxEnum);                          _assert_(vx_input);
	Input* vy_input             = element->GetInput(VyEnum);                          _assert_(vy_input);
	Input* vz_input             = element->GetInput(VzEnum);                          _assert_(vz_input);
	Input* enthalpy_input		 = element->GetInput(enthalpy_enum);					 _assert_(enthalpy_input);
	Input* pressure_input		 = element->GetInput(PressureEnum);							 _assert_(pressure_input);
	Input* watercolumn_input	 = element->GetInput(WatercolumnEnum);							 _assert_(watercolumn_input);
	Input* meltingrate_input	 = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);							 _assert_(meltingrate_input);
	Input* geothermalflux_input = element->GetInput(BasalforcingsGeothermalfluxEnum); _assert_(geothermalflux_input);
	IssmDouble  rho_ice			 = element->GetMaterialParameter(MaterialsRhoIceEnum);

	/*Build friction element, needed later: */
	Friction* friction=new Friction(element,3);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(4);
	Gauss* gaussup=element->NewGaussTop(4);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);
		gaussup->GaussPoint(ig);

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctions(basis,gauss);

		if(isdynamicbasalspc){
			enthalpy_input->GetInputValue(&enthalpy,gauss);
			enthalpy_input->GetInputValue(&enthalpyup,gaussup);
			pressure_input->GetInputValue(&pressure,gauss);
			pressure_input->GetInputValue(&pressureup,gaussup);
			watercolumn_input->GetInputValue(&watercolumn,gauss);
			meltingrate_input->GetInputValue(&meltingrate,gauss);
			state=GetThermalBasalCondition(element, enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate);
		}
		else
			state=0;

		switch (state) {
			case 0: case 1: case 2: case 3:
				// cold, dry base; cold, wet base; refreezing temperate base; thin temperate base: 
				// Apply basal surface forcing.
				// Interpolated values of enthalpy on gauss nodes may indicate cold base, 
				// although one node might have become temperate. So keep heat flux switched on.
				geothermalflux_input->GetInputValue(&geothermalflux,gauss);
				friction->GetAlpha2(&alpha2,gauss);
				vx_input->GetInputValue(&vx,gauss);
				vy_input->GetInputValue(&vy,gauss);
				vz_input->GetInputValue(&vz,gauss);
				basalfriction=alpha2*(vx*vx+vy*vy+vz*vz);
				heatflux=(basalfriction+geothermalflux)/(rho_ice);
				scalar=gauss->weight*Jdet*heatflux;
				if(dt!=0.) scalar=dt*scalar;
				for(i=0;i<numnodes;i++) 
					pe->values[i]+=scalar*basis[i];
				break;
			case 4:
				// temperate, thick melting base: set grad H*n=0
				for(i=0;i<numnodes;i++) 
					pe->values[i]+=0.;
				break;
			default:
				_printf0_("	unknown thermal basal state found!");
		}
	}

	/*Clean up and return*/
	delete gauss;
	delete gaussup;
	delete friction;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list_base);
	return pe;

}/*}}}*/
ElementVector* EnthalpyAnalysis::CreatePVectorShelf(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Get basal element*/
	if(!element->IsOnBase() || !element->IsFloating()) return NULL;

	IssmDouble  Hpmp,dt,Jdet,scalar_ocean,pressure;
	IssmDouble *xyz_list_base = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	Input*      pressure_input=element->GetInput(PressureEnum); _assert_(pressure_input);
	IssmDouble  gravity             = element->GetMaterialParameter(ConstantsGEnum);
	IssmDouble  rho_water           = element->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	IssmDouble  rho_ice             = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble  heatcapacity        = element->GetMaterialParameter(MaterialsHeatcapacityEnum);
	IssmDouble  mixed_layer_capacity= element->GetMaterialParameter(MaterialsMixedLayerCapacityEnum);
	IssmDouble  thermal_exchange_vel= element->GetMaterialParameter(MaterialsThermalExchangeVelocityEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(4);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctions(basis,gauss);

		pressure_input->GetInputValue(&pressure,gauss);
		Hpmp=element->PureIceEnthalpy(pressure);

		scalar_ocean=gauss->weight*Jdet*rho_water*mixed_layer_capacity*thermal_exchange_vel*Hpmp/(heatcapacity*rho_ice);
		if(reCast<bool,IssmDouble>(dt)) scalar_ocean=dt*scalar_ocean;

		for(int i=0;i<numnodes;i++) pe->values[i]+=scalar_ocean*basis[i];
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list_base);
	return pe;
}/*}}}*/
void           EnthalpyAnalysis::DrainWaterfraction(FemModel* femmodel){/*{{{*/
	/*Drain excess water fraction in ice column: */
	ComputeWaterfractionDrainage(femmodel);
	DrainageUpdateWatercolumn(femmodel);
	DrainageUpdateEnthalpy(femmodel);
}/*}}}*/
void				EnthalpyAnalysis::ComputeWaterfractionDrainage(FemModel* femmodel){/*{{{*/

	int i,k,numnodes;
	IssmDouble dt;
	Element* element= NULL;

	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	for(i=0;i<femmodel->elements->Size();i++){
		element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		numnodes=element->GetNumberOfNodes();
		IssmDouble* waterfractions= xNew<IssmDouble>(numnodes);
		IssmDouble* drainage= xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(waterfractions,WaterfractionEnum);
		for(k=0; k<numnodes;k++){
			drainage[k]=DrainageFunctionWaterfraction(waterfractions[k], dt);
		}
		element->AddInput(WaterfractionDrainageEnum,drainage,element->GetElementType());

		xDelete<IssmDouble>(waterfractions);
		xDelete<IssmDouble>(drainage);
	}
}/*}}}*/
void				EnthalpyAnalysis::DrainageUpdateWatercolumn(FemModel* femmodel){/*{{{*/

	int i,k,numnodes, numbasalnodes;
	IssmDouble dt;
	int* basalnodeindices=NULL;
	Element* element= NULL;
	
	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	/*depth-integrate the drained water fraction */
	femmodel->parameters->SetParam(WaterfractionDrainageEnum,InputToDepthaverageInEnum);
	femmodel->parameters->SetParam(WaterfractionDrainageIntegratedEnum,InputToDepthaverageOutEnum);
	depthaverage_core(femmodel);
	femmodel->parameters->SetParam(WaterfractionDrainageIntegratedEnum,InputToExtrudeEnum);
	extrudefrombase_core(femmodel);
	/*multiply depth-average by ice thickness*/
	for(i=0;i<femmodel->elements->Size();i++){
		element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		numnodes=element->GetNumberOfNodes();
		IssmDouble* drainage_int= xNew<IssmDouble>(numnodes);
		IssmDouble* thicknesses= xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(drainage_int,WaterfractionDrainageIntegratedEnum);
		element->GetInputListOnNodes(thicknesses,ThicknessEnum);
		for(k=0;k<numnodes;k++){
			drainage_int[k]*=thicknesses[k];
		}
		element->AddInput(WaterfractionDrainageIntegratedEnum, drainage_int, element->GetElementType());

		xDelete<IssmDouble>(drainage_int);
		xDelete<IssmDouble>(thicknesses);
	}

	/*update water column*/
	for(i=0;i<femmodel->elements->Size();i++){
		element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		/* Check if ice in element */
		if(!element->IsIceInElement()) continue;
		if(!element->IsOnBase()) continue; 

		numnodes=element->GetNumberOfNodes();
		IssmDouble* watercolumn= xNew<IssmDouble>(numnodes);
		IssmDouble* drainage_int= xNew<IssmDouble>(numnodes);
		element->GetInputListOnNodes(watercolumn,WatercolumnEnum);
		element->GetInputListOnNodes(drainage_int,WaterfractionDrainageIntegratedEnum);

		element->BasalNodeIndices(&numbasalnodes,&basalnodeindices,element->GetElementType());
		for(k=0;k<numbasalnodes;k++){
			watercolumn[basalnodeindices[k]]+=dt*drainage_int[basalnodeindices[k]];
		}
		element->AddInput(WatercolumnEnum, watercolumn, element->GetElementType());

		xDelete<IssmDouble>(watercolumn);
		xDelete<IssmDouble>(drainage_int);
	}
	xDelete<int>(basalnodeindices);
}/*}}}*/
void				EnthalpyAnalysis::DrainageUpdateEnthalpy(FemModel* femmodel){/*{{{*/

	int i,k,numnodes;
	IssmDouble dt;
	Element* element= NULL;
	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	for(i=0;i<femmodel->elements->Size();i++){
		element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		numnodes=element->GetNumberOfNodes();
		IssmDouble* enthalpies= xNew<IssmDouble>(numnodes);
		IssmDouble* pressures= xNew<IssmDouble>(numnodes);
		IssmDouble* temperatures= xNew<IssmDouble>(numnodes);
		IssmDouble* waterfractions= xNew<IssmDouble>(numnodes);
		IssmDouble* drainage= xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(pressures,PressureEnum);
		element->GetInputListOnNodes(temperatures,TemperatureEnum);
		element->GetInputListOnNodes(waterfractions,WaterfractionEnum);
		element->GetInputListOnNodes(drainage,WaterfractionDrainageEnum);

		for(k=0;k<numnodes;k++){
			waterfractions[k]-=dt*drainage[k];
			element->ThermalToEnthalpy(&enthalpies[k], temperatures[k], waterfractions[k], pressures[k]);
		}
		element->AddInput(WaterfractionEnum,waterfractions,element->GetElementType());
		element->AddInput(EnthalpyEnum,enthalpies,element->GetElementType());

		xDelete<IssmDouble>(enthalpies);
		xDelete<IssmDouble>(pressures);
		xDelete<IssmDouble>(temperatures);
		xDelete<IssmDouble>(waterfractions);
		xDelete<IssmDouble>(drainage);
	}
}/*}}}*/
IssmDouble     EnthalpyAnalysis::EnthalpyDiffusionParameter(Element* element,IssmDouble enthalpy,IssmDouble pressure){/*{{{*/

	IssmDouble heatcapacity             = element->GetMaterialParameter(MaterialsHeatcapacityEnum);
	IssmDouble temperateiceconductivity = element->GetMaterialParameter(MaterialsTemperateiceconductivityEnum);
	IssmDouble thermalconductivity      = element->GetMaterialParameter(MaterialsThermalconductivityEnum);

	if(enthalpy < PureIceEnthalpy(element,pressure)){
		return thermalconductivity/heatcapacity;
	}
	else{
		return temperateiceconductivity/heatcapacity;
	}
}/*}}}*/
IssmDouble     EnthalpyAnalysis::EnthalpyDiffusionParameterVolume(Element* element,int enthalpy_enum){/*{{{*/

	int         iv;
	IssmDouble  lambda;                   /* fraction of cold ice    */
	IssmDouble  kappa,kappa_c,kappa_t; /* enthalpy conductivities */
	IssmDouble  Hc,Ht;

	/*Get pressures and enthalpies on vertices*/
	int         numvertices = element->GetNumberOfVertices();
	IssmDouble* pressures   = xNew<IssmDouble>(numvertices);
	IssmDouble* enthalpies  = xNew<IssmDouble>(numvertices);
	IssmDouble* PIE         = xNew<IssmDouble>(numvertices);
	IssmDouble* dHpmp       = xNew<IssmDouble>(numvertices);
	element->GetInputListOnVertices(pressures,PressureEnum);
	element->GetInputListOnVertices(enthalpies,enthalpy_enum);
	for(iv=0;iv<numvertices;iv++){
		PIE[iv]   = PureIceEnthalpy(element,pressures[iv]);
		dHpmp[iv] = enthalpies[iv]-PIE[iv];
	}

	bool allequalsign = true;
	if(dHpmp[0]<0.){
		for(iv=1; iv<numvertices;iv++) allequalsign=(allequalsign && (dHpmp[iv]<0.));
	}
	else{
		for(iv=1; iv<numvertices;iv++) allequalsign=(allequalsign && (dHpmp[iv]>=0.));
	}

	if(allequalsign){
		kappa = EnthalpyDiffusionParameter(element,enthalpies[0],pressures[0]);
	}
	else{
		/* return harmonic mean of thermal conductivities, weighted by fraction of cold/temperate ice,
			cf Patankar 1980, pp44 */
		kappa_c = EnthalpyDiffusionParameter(element,PureIceEnthalpy(element,0.)-1.,0.);
		kappa_t = EnthalpyDiffusionParameter(element,PureIceEnthalpy(element,0.)+1.,0.);
		Hc=0.; Ht=0.;
		for(iv=0; iv<numvertices;iv++){
			if(enthalpies[iv]<PIE[iv])
			 Hc+=(PIE[iv]-enthalpies[iv]);
			else
			 Ht+=(enthalpies[iv]-PIE[iv]);
		}
		_assert_((Hc+Ht)>0.);
		lambda = Hc/(Hc+Ht);
		kappa  = kappa_c*kappa_t/(lambda*kappa_t+(1.-lambda)*kappa_c); // ==(lambda/kappa_c + (1.-lambda)/kappa_t)^-1
	}	

	/*Clean up and return*/
	xDelete<IssmDouble>(PIE);
	xDelete<IssmDouble>(dHpmp);
	xDelete<IssmDouble>(pressures);
	xDelete<IssmDouble>(enthalpies);
	return kappa;
}/*}}}*/
void           EnthalpyAnalysis::GetBAdvec(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3 B4 B5 B6] where Bi is of size 5*NDOF1. 
	 * For node i, Bi' can be expressed in the actual coordinate system
	 * by: 
	 *       Bi_advec =[ h ]
	 *                 [ h ]
	 *                 [ h ]
	 * where h is the interpolation function for node i.
	 *
	 * We assume B has been allocated already, of size: 3x(NDOF1*NUMNODESP1)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions*/
	IssmDouble* basis=xNew<IssmDouble>(numnodes);
	element->NodalFunctions(basis,gauss);

	/*Build B: */
	for(int i=0;i<numnodes;i++){
		B[numnodes*0+i] = basis[i];
		B[numnodes*1+i] = basis[i];
		B[numnodes*2+i] = basis[i];
	}

	/*Clean-up*/
	xDelete<IssmDouble>(basis);
}/*}}}*/
void           EnthalpyAnalysis::GetBAdvecprime(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3 B4 B5 B6] where Bi is of size 5*NDOF1. 
	 * For node i, Bi' can be expressed in the actual coordinate system
	 * by: 
	 *       Biprime_advec=[ dh/dx ]
	 *                     [ dh/dy ]
	 *                     [ dh/dz ]
	 * where h is the interpolation function for node i.
	 *
	 * We assume B has been allocated already, of size: 3x(NDOF1*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* dbasis=xNew<IssmDouble>(3*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B: */
	for(int i=0;i<numnodes;i++){
		B[numnodes*0+i] = dbasis[0*numnodes+i];
		B[numnodes*1+i] = dbasis[1*numnodes+i];
		B[numnodes*2+i] = dbasis[2*numnodes+i];
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);
}/*}}}*/
void           EnthalpyAnalysis::GetBasalConstraints(Vector<IssmDouble>* vec_spc,Element* element){/*{{{*/

	/*Intermediary*/
	bool        isdynamicbasalspc;
	IssmDouble	dt;

	/*Check wether dynamic basal boundary conditions are activated */
	element->FindParam(&isdynamicbasalspc,ThermalIsdynamicbasalspcEnum);
	if(!isdynamicbasalspc) return;

	element->FindParam(&dt,TimesteppingTimeStepEnum);
	if(dt==0.){
		GetBasalConstraintsSteadystate(vec_spc,element);
	}
	else{
		GetBasalConstraintsTransient(vec_spc,element);
	}
}/*}}}*/
void           EnthalpyAnalysis::GetBasalConstraintsSteadystate(Vector<IssmDouble>* vec_spc,Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return;

	/* Only update constraints at the base. 
	 * Floating ice is not affected by basal BC decision chart. */
	if(!(element->IsOnBase()) || element->IsFloating()) return;

	/*Intermediary*/
	int         numindices, numindicesup, state;
	int        *indices = NULL, *indicesup = NULL;
	IssmDouble	enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate;

	/*Get parameters and inputs: */
	Input* enthalpy_input		 = element->GetInput(EnthalpyPicardEnum);					 _assert_(enthalpy_input);
	Input* pressure_input		 = element->GetInput(PressureEnum);							 _assert_(pressure_input);
	Input* watercolumn_input	 = element->GetInput(WatercolumnEnum);							 _assert_(watercolumn_input);
	Input* meltingrate_input	 = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);							 _assert_(meltingrate_input);

	/*Fetch indices of basal & surface nodes for this finite element*/
	Penta *penta =  (Penta *) element; // TODO: add Basal-/SurfaceNodeIndices to element.h, and change this to Element*
	penta->BasalNodeIndices(&numindices,&indices,element->GetElementType());
	penta->SurfaceNodeIndices(&numindicesup,&indicesup,element->GetElementType());	_assert_(numindices==numindicesup);

	GaussPenta* gauss=new GaussPenta();
	GaussPenta* gaussup=new GaussPenta();
	for(int i=0;i<numindices;i++){
		gauss->GaussNode(element->GetElementType(),indices[i]);
		gaussup->GaussNode(element->GetElementType(),indicesup[i]);

		enthalpy_input->GetInputValue(&enthalpy,gauss);
		enthalpy_input->GetInputValue(&enthalpyup,gaussup);
		pressure_input->GetInputValue(&pressure,gauss);
		pressure_input->GetInputValue(&pressureup,gaussup);
		watercolumn_input->GetInputValue(&watercolumn,gauss);
		meltingrate_input->GetInputValue(&meltingrate,gauss);

		state=GetThermalBasalCondition(element, enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate);
		switch (state) {
			case 0:
				// cold, dry base: apply basal surface forcing
				vec_spc->SetValue(element->nodes[i]->Sid(),0.,INS_VAL);
				break;
			case 1:
				// cold, wet base: keep at pressure melting point 
				vec_spc->SetValue(element->nodes[i]->Sid(),1.,INS_VAL);
				break;
			case 2:
				// temperate, thin refreezing base: 
				vec_spc->SetValue(element->nodes[i]->Sid(),1.,INS_VAL);
				break;
			case 3:
				// temperate, thin melting base: set spc
				vec_spc->SetValue(element->nodes[i]->Sid(),1.,INS_VAL);
				break;
			case 4:
				// temperate, thick melting base:
				vec_spc->SetValue(element->nodes[i]->Sid(),1.,INS_VAL);
				break;
			default:
				_printf0_("	unknown thermal basal state found!");
		}
	}

	/*Free ressources:*/
	xDelete<int>(indices);
	xDelete<int>(indicesup);
	delete gauss;
	delete gaussup;
}/*}}}*/
void           EnthalpyAnalysis::GetBasalConstraintsTransient(Vector<IssmDouble>* vec_spc,Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return;

	/* Only update constraints at the base. 
	 * Floating ice is not affected by basal BC decision chart.*/
	if(!(element->IsOnBase()) || element->IsFloating()) return;

	/*Intermediary*/
	int         numindices, numindicesup, state;
	int        *indices = NULL, *indicesup = NULL;
	IssmDouble	enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate;

	/*Get parameters and inputs: */
	Input* enthalpy_input       = element->GetInput(EnthalpyEnum);                    _assert_(enthalpy_input); //TODO: check EnthalpyPicard?
	Input* pressure_input		 = element->GetInput(PressureEnum);							 _assert_(pressure_input);
	Input* watercolumn_input	 = element->GetInput(WatercolumnEnum);							 _assert_(watercolumn_input);
	Input* meltingrate_input	 = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);							 _assert_(meltingrate_input);

	/*Fetch indices of basal & surface nodes for this finite element*/
	Penta *penta =  (Penta *) element; // TODO: add Basal-/SurfaceNodeIndices to element.h, and change this to Element*
	penta->BasalNodeIndices(&numindices,&indices,element->GetElementType());
	penta->SurfaceNodeIndices(&numindicesup,&indicesup,element->GetElementType());	_assert_(numindices==numindicesup);

	GaussPenta* gauss=new GaussPenta();
	GaussPenta* gaussup=new GaussPenta();

	for(int i=0;i<numindices;i++){
		gauss->GaussNode(element->GetElementType(),indices[i]);
		gaussup->GaussNode(element->GetElementType(),indicesup[i]);
		
		enthalpy_input->GetInputValue(&enthalpy,gauss);
		enthalpy_input->GetInputValue(&enthalpyup,gaussup);
		pressure_input->GetInputValue(&pressure,gauss);
		pressure_input->GetInputValue(&pressureup,gaussup);
		watercolumn_input->GetInputValue(&watercolumn,gauss);
		meltingrate_input->GetInputValue(&meltingrate,gauss);

		state=GetThermalBasalCondition(element, enthalpy, enthalpyup, pressure, pressureup, watercolumn, meltingrate);

		switch (state) {
			case 0:
				// cold, dry base: apply basal surface forcing
				vec_spc->SetValue(element->nodes[i]->Sid(),0.,INS_VAL);
				break;
			case 1:
				// cold, wet base: keep at pressure melting point 
				vec_spc->SetValue(element->nodes[i]->Sid(),1.,INS_VAL);
				break;
			case 2:
				// temperate, thin refreezing base: release spc
				vec_spc->SetValue(element->nodes[i]->Sid(),0.,INS_VAL);
				break;
			case 3:
				// temperate, thin melting base: set spc
				vec_spc->SetValue(element->nodes[i]->Sid(),1.,INS_VAL);
				break;
			case 4:
				// temperate, thick melting base: set grad H*n=0
				vec_spc->SetValue(element->nodes[i]->Sid(),0.,INS_VAL);
				break;
			default:
				_printf0_("	unknown thermal basal state found!");
		}

	}

	/*Free ressources:*/
	xDelete<int>(indices);
	xDelete<int>(indicesup);
	delete gauss;
	delete gaussup;
}/*}}}*/
void           EnthalpyAnalysis::GetBConduct(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3 B4 B5 B6] where Bi is of size 5*NDOF1. 
	 * For node i, Bi' can be expressed in the actual coordinate system
	 * by: 
	 *       Bi_conduct=[ dh/dx ]
	 *                  [ dh/dy ]
	 *                  [ dh/dz ]
	 * where h is the interpolation function for node i.
	 *
	 * We assume B has been allocated already, of size: 3x(NDOF1*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* dbasis=xNew<IssmDouble>(3*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B: */
	for(int i=0;i<numnodes;i++){
		B[numnodes*0+i] = dbasis[0*numnodes+i];
		B[numnodes*1+i] = dbasis[1*numnodes+i];
		B[numnodes*2+i] = dbasis[2*numnodes+i];
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);
}/*}}}*/
void           EnthalpyAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,EnthalpyEnum);
}/*}}}*/
int            EnthalpyAnalysis::GetThermalBasalCondition(Element* element, IssmDouble enthalpy, IssmDouble enthalpyup, IssmDouble pressure, IssmDouble pressureup, IssmDouble watercolumn, IssmDouble meltingrate){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return -1;

	/* Only update Constraints at the base of grounded ice*/
	if(!(element->IsOnBase())) return -1;

	/*Intermediary*/
	int state=-1;
	IssmDouble	dt;

	/*Get parameters and inputs: */
	element->FindParam(&dt,TimesteppingTimeStepEnum);

	if(enthalpy<PureIceEnthalpy(element,pressure)){
		if(watercolumn<=0.) state=0; // cold, dry base
		else state=1; // cold, wet base (refreezing)
	}
	else{
		if(enthalpyup<PureIceEnthalpy(element,pressureup)){
			if((dt==0.) && (meltingrate<0.)) state=2;	// refreezing temperate base (non-physical, only for steadystate solver)
			else	state=3; // temperate base, but no temperate layer
		}
		else state=4; // temperate layer with positive thickness
	}

	_assert_(state>=0);
	return state;
}/*}}}*/
IssmDouble     EnthalpyAnalysis::GetWetIceConductivity(Element* element, IssmDouble enthalpy, IssmDouble pressure){/*{{{*/

	IssmDouble temperature, waterfraction;
	IssmDouble kappa_w = 0.6; // thermal conductivity of water (in W/m/K)
	IssmDouble kappa_i = element->GetMaterialParameter(MaterialsThermalconductivityEnum);
	element->EnthalpyToThermal(&temperature, &waterfraction, enthalpy, pressure);

	return (1.-waterfraction)*kappa_i + waterfraction*kappa_w;
}/*}}}*/
void           EnthalpyAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           EnthalpyAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	bool        converged;
	int         i,rheology_law;
	IssmDouble  B_average,s_average,T_average=0.,P_average=0.;
	IssmDouble  n=3.0;
	int        *doflist   = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = element->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	element->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values        = xNew<IssmDouble>(numnodes);
	IssmDouble* pressure      = xNew<IssmDouble>(numnodes);
	IssmDouble* surface       = xNew<IssmDouble>(numnodes);
	IssmDouble* B             = xNew<IssmDouble>(numnodes);
	IssmDouble* temperature   = xNew<IssmDouble>(numnodes);
	IssmDouble* waterfraction = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numnodes;i++){
		values[i]=solution[doflist[i]];

		/*Check solution*/
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*Get all inputs and parameters*/
	if(element->material->ObjectEnum()!=MatestarEnum) n=element->GetMaterialParameter(MaterialsRheologyNEnum);
	element->GetInputValue(&converged,ConvergedEnum);
	element->GetInputListOnNodes(&pressure[0],PressureEnum);
	if(converged){
		for(i=0;i<numnodes;i++){
			element->EnthalpyToThermal(&temperature[i],&waterfraction[i],values[i],pressure[i]);
			if(waterfraction[i]<0.) _error_("Negative water fraction found in solution vector");
			//if(waterfraction[i]>1.) _error_("Water fraction >1 found in solution vector");
		}
		element->AddInput(EnthalpyEnum,values,element->GetElementType());
		element->AddInput(WaterfractionEnum,waterfraction,element->GetElementType());
		element->AddInput(TemperatureEnum,temperature,element->GetElementType());

		/*Update Rheology only if converged (we must make sure that the temperature is below melting point
		 * otherwise the rheology could be negative*/
		element->FindParam(&rheology_law,MaterialsRheologyLawEnum);
		element->GetInputListOnNodes(&surface[0],SurfaceEnum);
		switch(rheology_law){
			case NoneEnum:
				/*Do nothing: B is not temperature dependent*/
				break;
			case BuddJackaEnum:
				for(i=0;i<numnodes;i++) B[i]=BuddJacka(temperature[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break;
			case CuffeyEnum:
				for(i=0;i<numnodes;i++) B[i]=Cuffey(temperature[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break;
			case CuffeyTemperateEnum:
				for(i=0;i<numnodes;i++) B[i]=CuffeyTemperate(temperature[i], waterfraction[i],n);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break;
			case PatersonEnum:
				for(i=0;i<numnodes;i++) B[i]=Paterson(temperature[i]);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break;
			case ArrheniusEnum:
				element->GetVerticesCoordinates(&xyz_list);
				for(i=0;i<numnodes;i++) B[i]=Arrhenius(temperature[i],surface[i]-xyz_list[i*3+2],n);
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType());
				break;
			case LliboutryDuvalEnum:
				for(i=0;i<numnodes;i++) B[i]=LliboutryDuval(values[i],pressure[i],n,element->GetMaterialParameter(MaterialsBetaEnum),element->GetMaterialParameter(ConstantsReferencetemperatureEnum),element->GetMaterialParameter(MaterialsHeatcapacityEnum),element->GetMaterialParameter(MaterialsLatentheatEnum)); 
				element->AddInput(MaterialsRheologyBEnum,&B[0],element->GetElementType()); 
				break; 
			default: _error_("Rheology law " << EnumToStringx(rheology_law) << " not supported yet");
		}
	}
	else{
		element->AddInput(EnthalpyPicardEnum,values,element->GetElementType());
	}

	/*Free ressources:*/
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(surface);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(temperature);
	xDelete<IssmDouble>(waterfraction);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
}/*}}}*/
void           EnthalpyAnalysis::PostProcessing(FemModel* femmodel){/*{{{*/

	/*Intermediaries*/
	bool computebasalmeltingrates=true;
	bool drainicecolumn=true;
	IssmDouble dt;

	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	if(drainicecolumn && (dt>0.))	DrainWaterfraction(femmodel);
	if(computebasalmeltingrates)	ComputeBasalMeltingrate(femmodel);

}/*}}}*/
IssmDouble     EnthalpyAnalysis::PureIceEnthalpy(Element* element,IssmDouble pressure){/*{{{*/

	IssmDouble heatcapacity         = element->GetMaterialParameter(MaterialsHeatcapacityEnum);
	IssmDouble referencetemperature = element->GetMaterialParameter(ConstantsReferencetemperatureEnum);

	return heatcapacity*(TMeltingPoint(element,pressure)-referencetemperature);
}/*}}}*/
IssmDouble     EnthalpyAnalysis::TMeltingPoint(Element* element,IssmDouble pressure){/*{{{*/

	IssmDouble meltingpoint = element->GetMaterialParameter(MaterialsMeltingpointEnum);
	IssmDouble beta         = element->GetMaterialParameter(MaterialsBetaEnum);

	return meltingpoint-beta*pressure;
}/*}}}*/
void           EnthalpyAnalysis::UpdateBasalConstraints(FemModel* femmodel){/*{{{*/

	/*Update basal dirichlet BCs for enthalpy: */
	Vector<IssmDouble>* spc           = NULL;
	IssmDouble*         serial_spc    = NULL;

	spc=new Vector<IssmDouble>(femmodel->nodes->NumberOfNodes(EnthalpyAnalysisEnum));
	/*First create a vector to figure out what elements should be constrained*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		GetBasalConstraints(spc,element);
	}

	/*Assemble and serialize*/
	spc->Assemble();
	serial_spc=spc->ToMPISerial();
	delete spc;

	/*Then update basal constraints nodes accordingly*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ApplyBasalConstraints(serial_spc,element);
	}

	femmodel->UpdateConstraintsx();

	/*Delete*/
	xDelete<IssmDouble>(serial_spc);
}/*}}}*/
void           EnthalpyAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
}/*}}}*/
