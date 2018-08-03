/*!\file Penta.cpp
 * \brief: implementation of the Penta object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*Element macros*/
#define NUMVERTICES   6
#define NUMVERTICES2D 3

/*Constructors/destructor/copy*/
Penta::~Penta(){/*{{{*/
	this->parameters=NULL;
}
/*}}}*/
Penta::Penta(int penta_id, int penta_sid, int index, IoModel* iomodel,int nummodels)/*{{{*/
	:ElementHook(nummodels,index+1,NUMVERTICES,iomodel){

	int penta_elements_ids[2];

	/*Checks in debugging mode*/
	_assert_(iomodel->Data("md.mesh.upperelements"));
	_assert_(iomodel->Data("md.mesh.lowerelements"));

	/*id: */
	this->id  = penta_id;
	this->sid = penta_sid;

	/*Build neighbors list*/
	if (xIsNan<IssmDouble>(iomodel->Data("md.mesh.upperelements")[index]) || iomodel->Data("md.mesh.upperelements")[index]==-1.) penta_elements_ids[1]=this->id; //upper penta is the same penta
	else                                    penta_elements_ids[1]=reCast<int,IssmDouble>((iomodel->Data("md.mesh.upperelements")[index]));
	if (xIsNan<IssmDouble>(iomodel->Data("md.mesh.lowerelements")[index]) || iomodel->Data("md.mesh.lowerelements")[index]==-1.) penta_elements_ids[0]=this->id; //lower penta is the same penta
	else                                    penta_elements_ids[0]=reCast<int,IssmDouble>((iomodel->Data("md.mesh.lowerelements")[index]));
	this->InitHookNeighbors(penta_elements_ids);

	//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
	this->parameters=NULL;

	/*intialize inputs: */
	this->inputs=new Inputs();

	/*initialize pointers:*/
	this->nodes             = NULL;
	this->vertices          = NULL;
	this->material          = NULL;
	this->matpar            = NULL;
	this->verticalneighbors = NULL;

	/*Only allocate pointer*/
	this->element_type_list=xNew<int>(nummodels);
}
/*}}}*/
Object* Penta::copy() {/*{{{*/

	int i;
	Penta* penta=NULL;

	penta=new Penta();

	//deal with PentaRef mother class
	int nanalyses = this->numanalyses;
	if(nanalyses > 0){
		penta->element_type_list=xNew<int>(nanalyses);
		for(i=0;i<nanalyses;i++) {
			if (this->element_type_list[i]) penta->element_type_list[i]=this->element_type_list[i];
			else penta->element_type_list[i] = 0;
		}
	}
	else penta->element_type_list = NULL;
	penta->element_type=this->element_type;
	penta->numanalyses=nanalyses;

	//deal with ElementHook mother class
	if (this->hnodes){
		penta->hnodes=xNew<Hook*>(penta->numanalyses);
		for(i=0;i<penta->numanalyses;i++){
			if (this->hnodes[i]) penta->hnodes[i] = (Hook*)(this->hnodes[i]->copy());
			else penta->hnodes[i] = NULL;
		}
	}
	else penta->hnodes = NULL;

	penta->hvertices = (Hook*)this->hvertices->copy();
	penta->hmaterial = (Hook*)this->hmaterial->copy();
	penta->hmatpar   = (Hook*)this->hmatpar->copy();
	if (this->hneighbors) penta->hneighbors = (Hook*)(this->hneighbors->copy());
	else penta->hneighbors = NULL;

	/*deal with Tria fields: */
	penta->id  = this->id;
	penta->sid = this->sid;
	if(this->inputs) penta->inputs = (Inputs*)(this->inputs->Copy());
	else penta->inputs=new Inputs();

	/*point parameters: */
	penta->parameters=this->parameters;

	/*recover objects: */
	if (this->nodes) {
		unsigned int num_nodes = 6;
		penta->nodes = xNew<Node*>(num_nodes); //we cannot rely on an analysis_counter to tell us which analysis_type we are running, so we just copy the nodes.
		for(i=0;i<num_nodes;i++) if(this->nodes[i]) penta->nodes[i]=this->nodes[i]; else penta->nodes[i] = NULL;
	}
	else penta->nodes = NULL;

	penta->vertices = (Vertex**)this->hvertices->deliverp();
	penta->material = (Material*)this->hmaterial->delivers();
	penta->matpar   = (Matpar*)this->hmatpar->delivers();
	penta->verticalneighbors = (Penta**)this->hneighbors->deliverp();

	return penta;

}
/*}}}*/
void Penta::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(PentaEnum);

	/*Call parent classes: */
	ElementHook::Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	Element::MarshallElement(pmarshalled_data,pmarshalled_data_size,marshall_direction,this->numanalyses);
	PentaRef::Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);

	vertices = (Vertex**)this->hvertices->deliverp();
	material = (Material*)this->hmaterial->delivers();
	matpar   = (Matpar*)this->hmatpar->delivers();
	verticalneighbors = (Penta**)this->hneighbors->deliverp();

}
/*}}}*/

/*Other*/
void       Penta::AddBasalInput(int input_enum,IssmDouble* values, int interpolation_enum){/*{{{*/

	_assert_(this->inputs);
	if(!IsOnBase()) return;
	else{
		if(interpolation_enum==P1Enum){
			int        i;
			IssmDouble extrudedvalues[NUMVERTICES];
			Penta*     penta=NULL;

			for(i=0;i<NUMVERTICES2D;i++){
				extrudedvalues[i]=values[i];
				extrudedvalues[i+NUMVERTICES2D]=values[i];
			}
			penta=this;
			for(;;){
				penta->inputs->AddInput(new PentaInput(input_enum,&extrudedvalues[0],P1Enum));
				if (penta->IsOnSurface()) break;
				penta=penta->GetUpperPenta(); _assert_(penta->Id()!=this->id);
			}
		}
		else _error_("not implemented yet");
	}
}
/*}}}*/
void       Penta::AddInput(int input_enum,IssmDouble* values, int interpolation_enum){/*{{{*/

	_assert_(this->inputs);
	this->inputs->AddInput(new PentaInput(input_enum,values,interpolation_enum));
}
/*}}}*/
void     Penta::BasalNodeIndices(int* pnumindices,int** pindices,int finiteelement){/*{{{*/

	PentaRef::BasalNodeIndices(pnumindices,pindices,finiteelement);

}
/*}}}*/
void       Penta::AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part){/*{{{*/
	_error_("Not supported yet!");
}
/*}}}*/
void       Penta::CalvingRateDev(){/*{{{*/

	if(!this->IsOnBase()) return;

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  calvingratex[NUMVERTICES];
	IssmDouble  calvingratey[NUMVERTICES];
	IssmDouble  calvingrate[NUMVERTICES];
	IssmDouble  lambda1,lambda2,ex,ey,vx,vy,vel;
	IssmDouble  sigma_vm,sigma_max,sigma_max_floating,sigma_max_grounded;
	IssmDouble  epse_2,groundedice;

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Depth average B for stress calculation*/
	this->InputDepthAverageAtBase(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Retrieve all inputs and parameters we will need*/
	Input* vx_input = inputs->GetInput(VxAverageEnum); _assert_(vx_input);
	Input* vy_input = inputs->GetInput(VyAverageEnum); _assert_(vy_input);
	Input* gr_input = inputs->GetInput(MaskGroundediceLevelsetEnum); _assert_(gr_input);
	IssmDouble  B   = this->GetMaterialParameter(MaterialsRheologyBbarEnum);
	IssmDouble  n   = this->GetMaterialParameter(MaterialsRheologyNEnum);
	this->parameters->FindParam(&sigma_max_floating,CalvingStressThresholdFloatingiceEnum);
	this->parameters->FindParam(&sigma_max_grounded,CalvingStressThresholdGroundediceEnum);

	/* Start looping on the number of vertices: */
	GaussPenta* gauss=new GaussPenta();
	for(int iv=0;iv<3;iv++){
		gauss->GaussVertex(iv);

		/*Get velocity components and thickness*/
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		gr_input->GetInputValue(&groundedice,gauss);
		vel=sqrt(vx*vx+vy*vy)+1.e-14;

		/*Compute strain rate and viscosity: */
		this->StrainRateSSA(&epsilon[0],&xyz_list[0][0],gauss,vx_input,vy_input);

		/*Get Eigen values*/
		Matrix2x2Eigen(&lambda1,&lambda2,&ex,&ey,epsilon[0],epsilon[2],epsilon[1]);
		_assert_(!xIsNan<IssmDouble>(lambda1));
		_assert_(!xIsNan<IssmDouble>(lambda2));

		/*Process Eigen values (only account for extension)*/
		lambda1 = max(lambda1,0.);
		lambda2 = max(lambda2,0.);

		/*Calculate sigma_vm*/
		epse_2    = 1./2. *(lambda1*lambda1 + lambda2*lambda2);
		sigma_vm  = sqrt(3.) * B * pow(epse_2,1./(2.*n));

		/*Tensile stress threshold*/
		if(groundedice<0)
		 sigma_max = sigma_max_floating;
		else
		 sigma_max = sigma_max_grounded;

		/*Assign values*/
		calvingratex[iv]=vx*sigma_vm/sigma_max;
		calvingratey[iv]=vy*sigma_vm/sigma_max;
		calvingrate[iv]=sqrt(calvingratex[iv]*calvingratex[iv] + calvingratey[iv]*calvingratey[iv]);
	}

	/*Add input*/
	this->inputs->AddInput(new PentaInput(CalvingratexEnum,&calvingratex[0],P1Enum));
	this->inputs->AddInput(new PentaInput(CalvingrateyEnum,&calvingratey[0],P1Enum));
	this->inputs->AddInput(new PentaInput(CalvingCalvingrateEnum,&calvingrate[0],P1Enum));

	this->InputExtrude(CalvingratexEnum,-1);
	this->InputExtrude(CalvingrateyEnum,-1);
	this->InputExtrude(CalvingCalvingrateEnum,-1);

	/*Clean up and return*/
	delete gauss;
}
/*}}}*/
void       Penta::CalvingRateLevermann(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	GaussPenta* gauss=NULL;
	IssmDouble  vx,vy,vel;
	IssmDouble  strainparallel;
	IssmDouble  propcoeff;
	IssmDouble  strainperpendicular;
	IssmDouble  calvingratex[NUMVERTICES];
	IssmDouble  calvingratey[NUMVERTICES];
	IssmDouble  calvingrate[NUMVERTICES];


	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs and parameters we will need*/
	Input* vx_input=inputs->GetInput(VxEnum);																		_assert_(vx_input);
	Input* vy_input=inputs->GetInput(VyEnum);																		_assert_(vy_input);
	Input* strainparallel_input=inputs->GetInput(StrainRateparallelEnum);								_assert_(strainparallel_input);
	Input* strainperpendicular_input=inputs->GetInput(StrainRateperpendicularEnum);              _assert_(strainperpendicular_input);
	Input* levermanncoeff_input=inputs->GetInput(CalvinglevermannCoeffEnum);                     _assert_(levermanncoeff_input);

	/* Start looping on the number of vertices: */
	gauss=new GaussPenta();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/* Get the value we need*/
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vel=vx*vx+vy*vy;
		strainparallel_input->GetInputValue(&strainparallel,gauss);
		strainperpendicular_input->GetInputValue(&strainperpendicular,gauss);
		levermanncoeff_input->GetInputValue(&propcoeff,gauss);

		/*Calving rate proportionnal to the positive product of the strain rate along the ice flow direction and the strain rate perpendicular to the ice flow */
		calvingrate[iv]=propcoeff*strainparallel*strainperpendicular;	
		if(calvingrate[iv]<0){
			calvingrate[iv]=0;
		}
		calvingratex[iv]=calvingrate[iv]*vx/(sqrt(vel)+1.e-14);
		calvingratey[iv]=calvingrate[iv]*vy/(sqrt(vel)+1.e-14);
	}

	/*Add input*/
	this->inputs->AddInput(new PentaInput(CalvingratexEnum,&calvingratex[0],P1Enum));
	this->inputs->AddInput(new PentaInput(CalvingrateyEnum,&calvingratey[0],P1Enum));
	this->inputs->AddInput(new PentaInput(CalvingCalvingrateEnum,&calvingrate[0],P1Enum));

	/*Clean up and return*/
	delete gauss;

}
/*}}}*/
void       Penta::ComputeBasalStress(Vector<IssmDouble>* sigma_b){/*{{{*/

	int         i,j;
	int         dofv[3]={0,1,2};
	int         dofp[1]={3};
	int         analysis_type,approximation;
	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  xyz_list_tria[3][3];
	IssmDouble  rho_ice,gravity,FSreconditioning;
	IssmDouble  pressure,viscosity,Jdet2d;
	IssmDouble  bed_normal[3];
	IssmDouble  basalforce[3] = {0.};
	IssmDouble  epsilon[6]; /* epsilon=[exx,eyy,ezz,exy,exz,eyz];*/
	IssmDouble  stresstensor[6]={0.0};
	IssmDouble  sigma_xx,sigma_yy,sigma_zz;
	IssmDouble  sigma_xy,sigma_xz,sigma_yz;
	IssmDouble  surface=0,value=0;
	GaussPenta* gauss;

	/*retrive parameters: */
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	inputs->GetInputValue(&approximation,ApproximationEnum);

	/*Check analysis_types*/
	if (analysis_type!=StressbalanceAnalysisEnum) _error_("Not supported yet!");
	if (approximation!=FSApproximationEnum) _error_("Not supported yet!");

	/*retrieve some parameters: */
	this->parameters->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);

	if(!IsOnBase()){
		//put zero
		sigma_b->SetValue(id-1,0.0,INS_VAL);
		return;
	}

	/*recovre material parameters: */
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	gravity=matpar->GetMaterialParameter(ConstantsGEnum);

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	for(i=0;i<3;i++) for(j=0;j<3;j++) xyz_list_tria[i][j]=xyz_list[i][j];

	/*Retrieve all inputs we will be needing: */
	Input* pressure_input=inputs->GetInput(PressureEnum); _assert_(pressure_input);
	Input* vx_input=inputs->GetInput(VxEnum);             _assert_(vx_input);
	Input* vy_input=inputs->GetInput(VyEnum);             _assert_(vy_input);
	Input* vz_input=inputs->GetInput(VzEnum);             _assert_(vz_input);

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussPenta(0,1,2,2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateFS(&epsilon[0],&xyz_list[0][0],gauss,vx_input,vy_input,vz_input);
		this->material->ViscosityFS(&viscosity,3,&xyz_list[0][0],gauss,vx_input,vy_input,vz_input);
		pressure_input->GetInputValue(&pressure,gauss);

		/*Compute Stress*/
		sigma_xx=2*viscosity*epsilon[0]-pressure*FSreconditioning; // sigma = nu eps - pressure
		sigma_yy=2*viscosity*epsilon[1]-pressure*FSreconditioning;
		sigma_zz=2*viscosity*epsilon[2]-pressure*FSreconditioning;
		sigma_xy=2*viscosity*epsilon[3];
		sigma_xz=2*viscosity*epsilon[4];
		sigma_yz=2*viscosity*epsilon[5];

		/*Get normal vector to the bed */
		NormalBase(&bed_normal[0],&xyz_list_tria[0][0]);

		/*basalforce*/
		basalforce[0] += sigma_xx*bed_normal[0] + sigma_xy*bed_normal[1] + sigma_xz*bed_normal[2];
		basalforce[1] += sigma_xy*bed_normal[0] + sigma_yy*bed_normal[1] + sigma_yz*bed_normal[2];
		basalforce[2] += sigma_xz*bed_normal[0] + sigma_yz*bed_normal[1] + sigma_zz*bed_normal[2];

		GetTriaJacobianDeterminant(&Jdet2d, &xyz_list_tria[0][0],gauss);
		value+=sigma_zz*Jdet2d*gauss->weight;
		surface+=Jdet2d*gauss->weight;
	}
	value=value/surface;

	/*Add value to output*/
	sigma_b->SetValue(id-1,value,INS_VAL);
}
/*}}}*/
void       Penta::ComputeDeviatoricStressTensor(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  viscosity;
	IssmDouble  epsilon[6]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  tau_xx[NUMVERTICES];
	IssmDouble	tau_yy[NUMVERTICES];
	IssmDouble	tau_zz[NUMVERTICES];
	IssmDouble  tau_xy[NUMVERTICES];
	IssmDouble	tau_xz[NUMVERTICES];
	IssmDouble	tau_yz[NUMVERTICES];
	IssmDouble	tau_eff[NUMVERTICES];
	GaussPenta* gauss=NULL;

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	Input* vx_input=inputs->GetInput(VxEnum);             _assert_(vx_input);
	Input* vy_input=inputs->GetInput(VyEnum);             _assert_(vy_input);
	Input* vz_input=inputs->GetInput(VzEnum);             _assert_(vz_input);

	/* Start looping on the number of vertices: */
	gauss=new GaussPenta();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateFS(&epsilon[0],&xyz_list[0][0],gauss,vx_input,vy_input,vz_input);
		this->material->ViscosityFS(&viscosity,3,&xyz_list[0][0],gauss,vx_input,vy_input,vz_input);

		/*Compute Stress*/
		tau_xx[iv]=2*viscosity*epsilon[0]; // tau = nu eps 
		tau_yy[iv]=2*viscosity*epsilon[1];
		tau_zz[iv]=2*viscosity*epsilon[2];
		tau_xy[iv]=2*viscosity*epsilon[3];
		tau_xz[iv]=2*viscosity*epsilon[4];
		tau_yz[iv]=2*viscosity*epsilon[5];

		tau_eff[iv] = tau_xx[iv]*tau_xx[iv] + tau_yy[iv]*tau_yy[iv] + tau_zz[iv]*tau_zz[iv] +
		  2*tau_xy[iv]*tau_xy[iv] + 2*tau_xz[iv]*tau_xz[iv] + 2*tau_yz[iv]*tau_yz[iv];

		tau_eff[iv] = sqrt(tau_eff[iv]/2.);
	}

	/*Add Stress tensor components into inputs*/
	this->inputs->AddInput(new PentaInput(DeviatoricStressxxEnum,&tau_xx[0],P1Enum));
	this->inputs->AddInput(new PentaInput(DeviatoricStressxyEnum,&tau_xy[0],P1Enum));
	this->inputs->AddInput(new PentaInput(DeviatoricStressxzEnum,&tau_xz[0],P1Enum));
	this->inputs->AddInput(new PentaInput(DeviatoricStressyyEnum,&tau_yy[0],P1Enum));
	this->inputs->AddInput(new PentaInput(DeviatoricStressyzEnum,&tau_yz[0],P1Enum));
	this->inputs->AddInput(new PentaInput(DeviatoricStresszzEnum,&tau_zz[0],P1Enum));
	this->inputs->AddInput(new PentaInput(DeviatoricStresseffectiveEnum,&tau_eff[0],P1Enum));

	/*Clean up and return*/
	delete gauss;
}
/*}}}*/
void       Penta::ComputeStressTensor(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  pressure,viscosity;
	IssmDouble  epsilon[6]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  sigma_xx[NUMVERTICES];
	IssmDouble	sigma_yy[NUMVERTICES];
	IssmDouble	sigma_zz[NUMVERTICES];
	IssmDouble  sigma_xy[NUMVERTICES];
	IssmDouble	sigma_xz[NUMVERTICES];
	IssmDouble	sigma_yz[NUMVERTICES];
	GaussPenta* gauss=NULL;

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	Input* pressure_input=inputs->GetInput(PressureEnum); _assert_(pressure_input);
	Input* vx_input=inputs->GetInput(VxEnum);             _assert_(vx_input);
	Input* vy_input=inputs->GetInput(VyEnum);             _assert_(vy_input);
	Input* vz_input=inputs->GetInput(VzEnum);             _assert_(vz_input);

	/* Start looping on the number of vertices: */
	gauss=new GaussPenta();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateFS(&epsilon[0],&xyz_list[0][0],gauss,vx_input,vy_input,vz_input);
		this->material->ViscosityFS(&viscosity,3,&xyz_list[0][0],gauss,vx_input,vy_input,vz_input);
		pressure_input->GetInputValue(&pressure,gauss);

		/*Compute Stress*/
		sigma_xx[iv]=2*viscosity*epsilon[0]-pressure; // sigma = nu eps - pressure
		sigma_yy[iv]=2*viscosity*epsilon[1]-pressure;
		sigma_zz[iv]=2*viscosity*epsilon[2]-pressure;
		sigma_xy[iv]=2*viscosity*epsilon[3];
		sigma_xz[iv]=2*viscosity*epsilon[4];
		sigma_yz[iv]=2*viscosity*epsilon[5];
	}

	/*Add Stress tensor components into inputs*/
	this->inputs->AddInput(new PentaInput(StressTensorxxEnum,&sigma_xx[0],P1Enum));
	this->inputs->AddInput(new PentaInput(StressTensorxyEnum,&sigma_xy[0],P1Enum));
	this->inputs->AddInput(new PentaInput(StressTensorxzEnum,&sigma_xz[0],P1Enum));
	this->inputs->AddInput(new PentaInput(StressTensoryyEnum,&sigma_yy[0],P1Enum));
	this->inputs->AddInput(new PentaInput(StressTensoryzEnum,&sigma_yz[0],P1Enum));
	this->inputs->AddInput(new PentaInput(StressTensorzzEnum,&sigma_zz[0],P1Enum));

	/*Clean up and return*/
	delete gauss;
}
/*}}}*/
void       Penta::Configure(Elements* elementsin, Loads* loadsin, Nodes* nodesin,Vertices* verticesin, Materials* materialsin, Parameters* parametersin){/*{{{*/

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
	this->hneighbors->configure(elementsin);

	/*Now, go pick up the objects inside the hooks: */
	if (this->hnodes[analysis_counter]) this->nodes=(Node**)this->hnodes[analysis_counter]->deliverp();
	else this->nodes=NULL;
	this->vertices          = (Vertex**)this->hvertices->deliverp();
	this->material          = (Material*)this->hmaterial->delivers();
	this->matpar            = (Matpar*)this->hmatpar->delivers();
	this->verticalneighbors = (Penta**)this->hneighbors->deliverp();

	/*point parameters to real dataset: */
	this->parameters=parametersin;

	/*get inputs configured too: */
	this->inputs->Configure(parameters);
}
/*}}}*/
void       Penta::ControlInputSetGradient(IssmDouble* gradient,int enum_type,int control_index){/*{{{*/

	int    vertexpidlist[NUMVERTICES];
	IssmDouble grad_list[NUMVERTICES];
	Input* grad_input=NULL;
	Input* input=NULL;

	if(enum_type==MaterialsRheologyBbarEnum){
		input=(Input*)inputs->GetInput(MaterialsRheologyBEnum);
	}
	else if(enum_type==DamageDbarEnum){
		input=(Input*)inputs->GetInput(DamageDEnum);
	}
	else{
		input=inputs->GetInput(enum_type);
	}
	if (!input) _error_("Input " << EnumToStringx(enum_type) << " not found");
	if (input->ObjectEnum()!=ControlInputEnum) _error_("Input " << EnumToStringx(enum_type) << " is not a ControlInput");

	GradientIndexing(&vertexpidlist[0],control_index);
	for(int i=0;i<NUMVERTICES;i++) grad_list[i]=gradient[vertexpidlist[i]];
	grad_input=new PentaInput(GradientEnum,grad_list,P1Enum);
	((ControlInput*)input)->SetGradient(grad_input);

}/*}}}*/
void       Penta::ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum){/*{{{*/

	Input* input=NULL;

	if(control_enum==MaterialsRheologyBbarEnum){
		input=(Input*)inputs->GetInput(MaterialsRheologyBEnum);
	}
	else if(control_enum==DamageDbarEnum){
		input=(Input*)inputs->GetInput(DamageDEnum);
	}
	else{
		input=inputs->GetInput(control_enum);
	}
	if (!input) _error_("Input " << EnumToStringx(control_enum) << " not found");
	if (input->ObjectEnum()!=ControlInputEnum) _error_("Input " << EnumToStringx(control_enum) << " is not a ControlInput");

	int         sidlist[NUMVERTICES];
	int         connectivity[NUMVERTICES];
	IssmPDouble values[NUMVERTICES];
	IssmPDouble gradients[NUMVERTICES]; 
	IssmDouble  value,gradient;

	this->GetVerticesConnectivityList(&connectivity[0]);
	this->GetVerticesSidList(&sidlist[0]);

	GaussPenta* gauss=new GaussPenta();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		((ControlInput*)input)->GetInputValue(&value,gauss);
		((ControlInput*)input)->GetGradientValue(&gradient,gauss);

		values[iv]    = reCast<IssmPDouble>(value)/reCast<IssmPDouble>(connectivity[iv]);
		gradients[iv] = reCast<IssmPDouble>(gradient)/reCast<IssmPDouble>(connectivity[iv]);
	}
	delete gauss;

	vector_control->SetValues(NUMVERTICES,&sidlist[0],&values[0],ADD_VAL);
	vector_gradient->SetValues(NUMVERTICES,&sidlist[0],&gradients[0],ADD_VAL);

}/*}}}*/
void       Penta::ElementResponse(IssmDouble* presponse,int response_enum){/*{{{*/

	switch(response_enum){
		case MaterialsRheologyBbarEnum:
			*presponse=this->material->GetBbar();
			break;
		case DamageDbarEnum:
			*presponse=this->material->GetDbar();
			break;
		case VelEnum:
			{

				/*Get input:*/
				IssmDouble vel;
				Input* vel_input;

				vel_input=this->inputs->GetInput(VelEnum); _assert_(vel_input);
				vel_input->GetInputAverage(&vel);

				/*Assign output pointers:*/
				*presponse=vel;
			}
			break;
		default:  
			_error_("Response type " << EnumToStringx(response_enum) << " not supported yet!");
	}

}
/*}}}*/
void       Penta::ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz){/*{{{*/

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
int        Penta::FiniteElement(void){/*{{{*/
	return this->element_type;
}
/*}}}*/
IssmDouble Penta::FloatingArea(void){/*{{{*/

	/*Intermediaries*/
	int         domaintype;
	IssmDouble  phi,base_area;
	IssmDouble  xyz_list[NUMVERTICES][3];

	if(!IsIceInElement() || !IsOnBase())return 0.;

	/*Get problem dimension*/
	this->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain3DEnum) _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	phi=this->GetGroundedPortion(&xyz_list[0][0]);
	base_area= 1./2.*fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	/*Clean up and return*/
	return (1-phi)*base_area;
}
/*}}}*/
void       Penta::FSContactMigration(Vector<IssmDouble>* vertexgrounded,Vector<IssmDouble>* vertexfloating){/*{{{*/

	if(!IsOnBase()) return;

	int approximation;
	inputs->GetInputValue(&approximation,ApproximationEnum);
	if(approximation==HOApproximationEnum || approximation==SSAApproximationEnum || approximation==SSAHOApproximationEnum){
		for(int i=0;i<NUMVERTICES;i++){
			vertexgrounded->SetValue(vertices[i]->Pid(),+9999.,INS_VAL);
			vertexfloating->SetValue(vertices[i]->Pid(),+9999.,INS_VAL);
		}
	}
	else {
		/*Intermediaries*/
		IssmDouble* xyz_list = NULL;
		IssmDouble  pressure,water_pressure,sigma_nn,viscosity,bed,base;
		IssmDouble  bed_normal[3];
		IssmDouble  epsilon[6]; /* epsilon=[exx eyy ezz exy exz eyz];*/
		IssmDouble  surface=0,value=0;
		bool grounded;

		/* Get node coordinates and dof list: */
		GetVerticesCoordinates(&xyz_list);

		/*Retrieve all inputs we will be needing: */
		Input* pressure_input = inputs->GetInput(PressureEnum); _assert_(pressure_input);
		Input* base_input     = inputs->GetInput(BaseEnum);     _assert_(base_input);
		Input* bed_input      = inputs->GetInput(BedEnum);      _assert_(bed_input);
		Input* vx_input       = inputs->GetInput(VxEnum);       _assert_(vx_input);
		Input* vy_input       = inputs->GetInput(VyEnum);       _assert_(vy_input);
		Input* vz_input       = inputs->GetInput(VzEnum);       _assert_(vz_input);

		/*Create gauss point in the middle of the basal edge*/
		Gauss* gauss=NewGaussBase(1);
		gauss->GaussPoint(0);

		if(!IsFloating()){ 
			/*Check for basal force only if grounded and touching GL*/
			this->StrainRateFS(&epsilon[0],xyz_list,gauss,vx_input,vy_input,vz_input);
			this->material->ViscosityFS(&viscosity,3,xyz_list,gauss,vx_input,vy_input,vz_input);
			pressure_input->GetInputValue(&pressure, gauss);
			base_input->GetInputValue(&base, gauss); _assert_(base<0.);

			/*Compute Stress*/
			IssmDouble sigma_xx=2.*viscosity*epsilon[0]-pressure;
			IssmDouble sigma_yy=2.*viscosity*epsilon[1]-pressure;
			IssmDouble sigma_zz=2.*viscosity*epsilon[2]-pressure;
			IssmDouble sigma_xy=2.*viscosity*epsilon[3];
			IssmDouble sigma_xz=2.*viscosity*epsilon[4];
			IssmDouble sigma_yz=2.*viscosity*epsilon[5];

			/*Get normal vector to the bed */
			NormalBase(&bed_normal[0],xyz_list);

			/*basalforce*/
			sigma_nn = sigma_xx*bed_normal[0]*bed_normal[0] + sigma_yy*bed_normal[1]*bed_normal[1] + sigma_zz*bed_normal[2]*bed_normal[2] 
			  + 2.*sigma_xy*bed_normal[0]*bed_normal[1] + 2.*sigma_xz*bed_normal[0]*bed_normal[2] + 2.*sigma_yz*bed_normal[1]*bed_normal[2];

			/*Compute water pressure*/
			IssmDouble rho_ice   = matpar->GetMaterialParameter(MaterialsRhoIceEnum);
			IssmDouble rho_water = matpar->GetMaterialParameter(MaterialsRhoSeawaterEnum);
			IssmDouble gravity   = matpar->GetMaterialParameter(ConstantsGEnum);
			water_pressure=gravity*rho_water*base;

			/*Compare basal stress to water pressure and determine whether it should ground*/
			if (sigma_nn<water_pressure) grounded=true;
			else                         grounded=false;
		}
		else{
			/*Check for basal elevation if floating*/
			base_input->GetInputValue(&base, gauss);
			bed_input->GetInputValue(&bed, gauss);
			if(base<bed) grounded=true;
			else          grounded=false;
		}
		for(int i=0;i<NUMVERTICES;i++){
			if(grounded) vertexgrounded->SetValue(vertices[i]->Pid(),+1.,INS_VAL);
			else         vertexfloating->SetValue(vertices[i]->Pid(),+1.,INS_VAL);
		}

		/*clean up*/
		delete gauss;
		xDelete<IssmDouble>(xyz_list);
	}
}
/*}}}*/
void       Penta::GetAreaCoordinates(IssmDouble* area_coordinates,IssmDouble* xyz_zero,IssmDouble* xyz_list,int numpoints){/*{{{*/
	/*Computeportion of the element that is grounded*/ 

	int         i,j,k;
	IssmDouble  area_init,area_portion;
	IssmDouble  xyz_bis[3][3];

	area_init=fabs(xyz_list[1*3+0]*xyz_list[2*3+1] - xyz_list[1*3+1]*xyz_list[2*3+0] + xyz_list[0*3+0]*xyz_list[1*3+1] - xyz_list[0*3+1]*xyz_list[1*3+0] + xyz_list[2*3+0]*xyz_list[0*3+1] - xyz_list[2*3+1]*xyz_list[0*3+0])/2.;

	/*Initialize xyz_list with original xyz_list of triangle coordinates*/
	for(j=0;j<3;j++){ 
		for(k=0;k<3;k++){
			xyz_bis[j][k]=xyz_list[j*3+k];
		}
	}
	for(i=0;i<numpoints;i++){
		for(j=0;j<3;j++){ 
			for(k=0;k<3;k++){
				/*Change appropriate line*/
				xyz_bis[j][k]=xyz_zero[i*3+k];
			}

			/*Compute area fraction*/
			area_portion=fabs(xyz_bis[1][0]*xyz_bis[2][1] - xyz_bis[1][1]*xyz_bis[2][0] + xyz_bis[0][0]*xyz_bis[1][1] - xyz_bis[0][1]*xyz_bis[1][0] + xyz_bis[2][0]*xyz_bis[0][1] - xyz_bis[2][1]*xyz_bis[0][0])/2.;
			*(area_coordinates+3*i+j)=area_portion/area_init;

			/*Reinitialize xyz_list*/
			for(k=0;k<3;k++){
				/*Reinitialize xyz_list with original coordinates*/
				xyz_bis[j][k]=xyz_list[j*3+k];
			}
		}
	}
}
/*}}}*/
Element*   Penta::GetBasalElement(void){/*{{{*/

	/*Output*/
	Element* element=this->GetBasalPenta();
	return element;
}
/*}}}*/
Penta*     Penta::GetBasalPenta(void){/*{{{*/

	/*Output*/
	Penta* penta=NULL;

	/*Go through all pentas till the bed is reached*/
	penta=this;
	for(;;){
		/*Stop if we have reached the surface, else, take lower penta*/
		if (penta->IsOnBase()) break;

		/* get lower Penta*/
		penta=penta->GetLowerPenta();
		_assert_(penta->Id()!=this->id);
	}

	/*return output*/
	return penta;
}
/*}}}*/
int        Penta::GetElementType(){/*{{{*/

	/*return PentaRef field*/
	return this->element_type;
}
/*}}}*/
void       Penta::GetGroundedPart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlyfloating){/*{{{*/
	/*Computeportion of the element that is grounded*/ 

	bool               floating=true;
	int                point;
	const IssmPDouble  epsilon= 1.e-15;
	IssmDouble         gl[NUMVERTICES];
	IssmDouble         f1,f2;

	/*Recover parameters and values*/
	GetInputListOnVertices(&gl[0],MaskGroundediceLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	/*Check that not all nodes are grounded or floating*/
	if(gl[0]>0 && gl[1]>0 && gl[2]>0){ // All grounded
		point=0;
		f1=1.;
		f2=1.;
	}
	else if(gl[0]<0 && gl[1]<0 && gl[2]<0){ //All floating
		point=0;
		f1=0.;
		f2=0.;
	}
	else{
		if(gl[0]*gl[1]*gl[2]<0) floating=false;

		if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2
			point=2;
			f1=gl[2]/(gl[2]-gl[0]);
			f2=gl[2]/(gl[2]-gl[1]);
		}
		else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2
			point=0;
			f1=gl[0]/(gl[0]-gl[1]);
			f2=gl[0]/(gl[0]-gl[2]);
		}
		else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2
			point=1;
			f1=gl[1]/(gl[1]-gl[2]);
			f2=gl[1]/(gl[1]-gl[0]);
		}
	}
	*point1=point;
	*fraction1=f1;
	*fraction2=f2;
	*mainlyfloating=floating;
}
/*}}}*/
IssmDouble Penta::GetGroundedPortion(IssmDouble* xyz_list){/*{{{*/
	/*Computeportion of the element that is grounded*/ 

	bool               mainlyfloating = true;
	const IssmPDouble  epsilon= 1.e-15;
	IssmDouble         phi,s1,s2,area_init,area_grounded;
	IssmDouble         gl[NUMVERTICES];
	IssmDouble         xyz_bis[NUMVERTICES2D][3];

	/*Recover parameters and values*/
	GetInputListOnVertices(&gl[0],MaskGroundediceLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	/*Check that not all nodes are grounded or floating*/
	if(gl[0]>0 && gl[1]>0 && gl[2]>0){ // All grounded
		phi=1;
	}
	else if(gl[0]<0 && gl[1]<0 && gl[2]<0){ //All floating
		phi=0;
	}
	else{
		/*Figure out if two nodes are floating or grounded*/
		if(gl[0]*gl[1]*gl[2]>0) mainlyfloating=false;

		if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2
			/*Coordinates of point 2: same as initial point 2*/
			xyz_bis[2][0]=xyz_list[3*2+0];
			xyz_bis[2][1]=xyz_list[3*2+1];
			xyz_bis[2][2]=xyz_list[3*2+2];

			/*Portion of the segments*/
			s1=gl[2]/(gl[2]-gl[1]);
			s2=gl[2]/(gl[2]-gl[0]);

			/*New point 1*/
			xyz_bis[1][0]=xyz_list[3*2+0]+s1*(xyz_list[3*1+0]-xyz_list[3*2+0]);
			xyz_bis[1][1]=xyz_list[3*2+1]+s1*(xyz_list[3*1+1]-xyz_list[3*2+1]);
			xyz_bis[1][2]=xyz_list[3*2+2]+s1*(xyz_list[3*1+2]-xyz_list[3*2+2]);

			/*New point 0*/
			xyz_bis[0][0]=xyz_list[3*2+0]+s2*(xyz_list[3*0+0]-xyz_list[3*2+0]);
			xyz_bis[0][1]=xyz_list[3*2+1]+s2*(xyz_list[3*0+1]-xyz_list[3*2+1]);
			xyz_bis[0][2]=xyz_list[3*2+2]+s2*(xyz_list[3*0+2]-xyz_list[3*2+2]);
		}
		else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2
			/*Coordinates of point 0: same as initial point 2*/
			xyz_bis[0][0]=*(xyz_list+3*0+0);
			xyz_bis[0][1]=*(xyz_list+3*0+1);
			xyz_bis[0][2]=*(xyz_list+3*0+2);

			/*Portion of the segments*/
			s1=gl[0]/(gl[0]-gl[1]);
			s2=gl[0]/(gl[0]-gl[2]);

			/*New point 1*/
			xyz_bis[1][0]=*(xyz_list+3*0+0)+s1*(*(xyz_list+3*1+0)-*(xyz_list+3*0+0));
			xyz_bis[1][1]=*(xyz_list+3*0+1)+s1*(*(xyz_list+3*1+1)-*(xyz_list+3*0+1));
			xyz_bis[1][2]=*(xyz_list+3*0+2)+s1*(*(xyz_list+3*1+2)-*(xyz_list+3*0+2));

			/*New point 2*/
			xyz_bis[2][0]=*(xyz_list+3*0+0)+s2*(*(xyz_list+3*2+0)-*(xyz_list+3*0+0));
			xyz_bis[2][1]=*(xyz_list+3*0+1)+s2*(*(xyz_list+3*2+1)-*(xyz_list+3*0+1));
			xyz_bis[2][2]=*(xyz_list+3*0+2)+s2*(*(xyz_list+3*2+2)-*(xyz_list+3*0+2));
		}
		else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2
			/*Coordinates of point 1: same as initial point 2*/
			xyz_bis[1][0]=*(xyz_list+3*1+0);
			xyz_bis[1][1]=*(xyz_list+3*1+1);
			xyz_bis[1][2]=*(xyz_list+3*1+2);

			/*Portion of the segments*/
			s1=gl[1]/(gl[1]-gl[0]);
			s2=gl[1]/(gl[1]-gl[2]);

			/*New point 0*/
			xyz_bis[0][0]=*(xyz_list+3*1+0)+s1*(*(xyz_list+3*0+0)-*(xyz_list+3*1+0));
			xyz_bis[0][1]=*(xyz_list+3*1+1)+s1*(*(xyz_list+3*0+1)-*(xyz_list+3*1+1));
			xyz_bis[0][2]=*(xyz_list+3*1+2)+s1*(*(xyz_list+3*0+2)-*(xyz_list+3*1+2));

			/*New point 2*/
			xyz_bis[2][0]=*(xyz_list+3*1+0)+s2*(*(xyz_list+3*2+0)-*(xyz_list+3*1+0));
			xyz_bis[2][1]=*(xyz_list+3*1+1)+s2*(*(xyz_list+3*2+1)-*(xyz_list+3*1+1));
			xyz_bis[2][2]=*(xyz_list+3*1+2)+s2*(*(xyz_list+3*2+2)-*(xyz_list+3*1+2));
		}

		/*Compute fraction of grounded element*/
		GetTriaJacobianDeterminant(&area_init, xyz_list,NULL);
		GetTriaJacobianDeterminant(&area_grounded, &xyz_bis[0][0],NULL);
		if(mainlyfloating==true) area_grounded=area_init-area_grounded;
		phi=area_grounded/area_init;
	}

	if(phi>1. || phi<0.) _error_("Error. Problem with portion of grounded element: value should be between 0 and 1");

	return phi;
}
/*}}}*/
void       Penta::GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum){/*{{{*/
	
	/* Intermediaries */
	const int dim=3;
	int i, dir,nrfrontnodes;
	IssmDouble  levelset[NUMVERTICES];

	/*Recover parameters and values*/
	GetInputListOnVertices(&levelset[0],levelsetenum);

	int* indicesfront = xNew<int>(NUMVERTICES);
	/* Get basal nodes where there is no ice */
	nrfrontnodes=0;
	for(i=0;i<NUMVERTICES2D;i++){
		if(levelset[i]>=0.){
			indicesfront[nrfrontnodes]=i;
			nrfrontnodes++;
		}
	}
	_assert_(nrfrontnodes==2);

	/* arrange order of basal frontnodes such that they are oriented counterclockwise */
	if((NUMVERTICES2D+indicesfront[0]-indicesfront[1])%NUMVERTICES2D!=NUMVERTICES2D-1){
		int index=indicesfront[0];
		indicesfront[0]=indicesfront[1];
		indicesfront[1]=index;
	}	

	IssmDouble* xyz_front = xNew<IssmDouble>(2*dim*nrfrontnodes);
	/* Return basal and top front nodes */
	for(i=0;i<nrfrontnodes;i++){
		for(dir=0;dir<dim;dir++){
			int ind1=i*dim+dir, ind2=(2*nrfrontnodes-1-i)*dim+dir; // vertex structure front segment: base0, base1, top1, top0
			xyz_front[ind1]=xyz_list[dim*indicesfront[i]+dir];
			xyz_front[ind2]=xyz_list[dim*(indicesfront[i]+NUMVERTICES2D)+dir];
		}
	}

	*pxyz_front=xyz_front;

	xDelete<int>(indicesfront);
}/*}}}*/
void       Penta::GetInputValue(IssmDouble* pvalue,Node* node,int enumtype){/*{{{*/

	Input* input=inputs->GetInput(enumtype);
	if(!input) _error_("No input of type " << EnumToStringx(enumtype) << " found in tria");

	GaussPenta* gauss=new GaussPenta();
	gauss->GaussVertex(this->GetNodeIndex(node));

	input->GetInputValue(pvalue,gauss);
	delete gauss;
}
/*}}}*/
Penta*     Penta::GetLowerPenta(void){/*{{{*/

	Penta* lower_penta=NULL;

	lower_penta=(Penta*)verticalneighbors[0]; //first one (0) under, second one (1) above

	return lower_penta;
}
/*}}}*/
Node*      Penta::GetNode(int node_number){/*{{{*/
	_assert_(node_number>=0); 
	_assert_(node_number<this->NumberofNodes(this->element_type)); 
	return this->nodes[node_number];
}
/*}}}*/
int        Penta::GetNodeIndex(Node* node){/*{{{*/

	_assert_(nodes);
	int numnodes = this->NumberofNodes(this->element_type);

	for(int i=0;i<numnodes;i++){
		if(node==nodes[i]) return i;
	}
	_error_("Node provided not found among element nodes");

}
/*}}}*/
int        Penta::GetNumberOfNodes(void){/*{{{*/
	return this->NumberofNodes(this->element_type);
}
/*}}}*/
int        Penta::GetNumberOfNodes(int enum_type){/*{{{*/
	return this->NumberofNodes(enum_type);
}
/*}}}*/
int        Penta::GetNumberOfVertices(void){/*{{{*/
	return NUMVERTICES; 
}
/*}}}*/
void       Penta::GetSolutionFromInputsOneDof(Vector<IssmDouble>* solution, int enum_type){/*{{{*/

	const int    numdof=NDOF1*NUMVERTICES;

	int          i;
	int*         doflist=NULL;
	IssmDouble   values[numdof];
	IssmDouble   enum_value;
	GaussPenta   *gauss=NULL;

	/*Get dof list: */
	GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	Input* enum_input=inputs->GetInput(enum_type); _assert_(enum_input);

	gauss=new GaussPenta();
	for(i=0;i<NUMVERTICES;i++){
		/*Recover temperature*/
		gauss->GaussVertex(i);
		enum_input->GetInputValue(&enum_value,gauss);
		values[i]=enum_value;
	}

	/*Add value to global vector*/
	solution->SetValues(numdof,doflist,values,INS_VAL);

	/*Free ressources:*/
	delete gauss;
	xDelete<int>(doflist);
}
/*}}}*/
Penta*     Penta::GetSurfacePenta(void){/*{{{*/

	/*Output*/
	Penta* penta=NULL;

	/*Go through all pentas till the surface is reached*/
	penta=this;
	for(;;){
		/*Stop if we have reached the surface, else, take upper penta*/
		if (penta->IsOnSurface()) break;

		/* get upper Penta*/
		penta=penta->GetUpperPenta();
		_assert_(penta->Id()!=this->id);
	}

	/*return output*/
	return penta;
}
/*}}}*/
Element*   Penta::GetUpperElement(void){/*{{{*/

	/*Output*/
	Element* upper_element=this->GetUpperPenta();
	return upper_element;
}
/*}}}*/
Penta*     Penta::GetUpperPenta(void){/*{{{*/

	Penta* upper_penta=NULL;

	upper_penta=(Penta*)verticalneighbors[1]; //first one (0) under, second one (1) above

	return upper_penta;
}
/*}}}*/
void       Penta::GetVectorFromControlInputs(Vector<IssmDouble>* vector,int control_enum,int control_index,const char* data,bool onsid){/*{{{*/

	int vertexidlist[NUMVERTICES];

	/*Get out if this is not an element input*/
	if(!IsInput(control_enum)) _error_("Enum "<<EnumToStringx(control_enum)<<" is not in IsInput");

	/*Prepare index list*/
	GradientIndexing(&vertexidlist[0],control_index,onsid);

	/*Get input (either in element or material)*/
	if(control_enum==MaterialsRheologyBbarEnum) control_enum=MaterialsRheologyBEnum;
	Input* input=inputs->GetInput(control_enum);
	if(!input) _error_("Input " << EnumToStringx(control_enum) << " not found in element");

	/*Check that it is a ControlInput*/
	if (input->ObjectEnum()!=ControlInputEnum){
		_error_("input " << EnumToStringx(control_enum) << " is not a ControlInput");
	}

	((ControlInput*)input)->GetVectorFromInputs(vector,&vertexidlist[0],data);
}
/*}}}*/
void       Penta::GetVerticesCoordinatesBase(IssmDouble** pxyz_list){/*{{{*/

	IssmDouble* xyz_list = xNew<IssmDouble>(NUMVERTICES2D*3);
	::GetVerticesCoordinates(xyz_list,this->vertices,NUMVERTICES2D);

	/*Assign output pointer*/
	*pxyz_list = xyz_list;

}/*}}}*/
void       Penta::GetVerticesCoordinatesTop(IssmDouble** pxyz_list){/*{{{*/

	IssmDouble* xyz_list = xNew<IssmDouble>(NUMVERTICES2D*3);
	::GetVerticesCoordinates(xyz_list,&this->vertices[3],NUMVERTICES2D);

	/*Assign output pointer*/
	*pxyz_list = xyz_list;

}/*}}}*/
IssmDouble Penta::GroundedArea(void){/*{{{*/

	/*Intermediaries*/
	int         domaintype;
	IssmDouble  phi,base_area;
	IssmDouble  xyz_list[NUMVERTICES][3];

	if(!IsIceInElement() || !IsOnBase())return 0.;

	/*Get problem dimension*/
	this->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain3DEnum) _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	phi=this->GetGroundedPortion(&xyz_list[0][0]);
	base_area= 1./2.*fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	/*Clean up and return*/
	return phi*base_area;
}
/*}}}*/
IssmDouble Penta::IceMass(void){/*{{{*/

	IssmDouble rho_ice; 
	
	if(!IsIceInElement())return 0.; //do not contribute to the volume of the ice!

	/*recover ice density: */
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);

	return rho_ice*this->IceVolume();
}
/*}}}*/
IssmDouble Penta::IceVolume(void){/*{{{*/

	/*The volume of a troncated prism is base * 1/3 sum(length of edges)*/
	IssmDouble base,height;
	IssmDouble xyz_list[NUMVERTICES][3];

	if(!IsIceInElement())return 0;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Pentangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2.*fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	/*Now get the average height*/
	height = 1./3.*((xyz_list[3][2]-xyz_list[0][2])+(xyz_list[4][2]-xyz_list[1][2])+(xyz_list[5][2]-xyz_list[2][2]));

	/*Return: */
	return base*height;
}
/*}}}*/
IssmDouble Penta::IceVolumeAboveFloatation(void){/*{{{*/

	/*Volume above floatation: H + rho_water/rho_ice*bathymetry for nodes on the bed*/
	IssmDouble rho_ice,rho_water;
	IssmDouble base,bed,surface,bathymetry;
	IssmDouble xyz_list[NUMVERTICES][3];

	if(!IsIceInElement() || IsFloating() || !IsOnBase())return 0;

	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	rho_water=matpar->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Pentangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2.*fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	/*Now get the average height above floatation*/
	Input* surface_input    = inputs->GetInput(SurfaceEnum);    _assert_(surface_input);
	Input* base_input        = inputs->GetInput(BaseEnum);        _assert_(base_input);
	Input* bed_input = inputs->GetInput(BedEnum); _assert_(bed_input);
	surface_input->GetInputAverage(&surface);
	base_input->GetInputAverage(&bed);
	bed_input->GetInputAverage(&bathymetry);

	/*Return: */
	return base*(surface - bed + min( rho_water/rho_ice * bathymetry, 0.) );
}
/*}}}*/
void       Penta::InputControlUpdate(IssmDouble scalar,bool save_parameter){/*{{{*/

	/*Intermediary*/
	int    num_controls;
	int*   control_type=NULL;
	Input* input=NULL;

	/*retrieve some parameters: */
	this->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	this->parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);

	for(int i=0;i<num_controls;i++){

		if(control_type[i]==MaterialsRheologyBbarEnum){
			if (!IsOnBase()) goto cleanup_and_return;
			input=(Input*)this->inputs->GetInput(MaterialsRheologyBEnum); _assert_(input);
		}
		else if(control_type[i]==DamageDbarEnum){
			if (!IsOnBase()) goto cleanup_and_return;
			input=(Input*)this->inputs->GetInput(DamageDEnum); _assert_(input);
		}
		else{
			input=(Input*)this->inputs->GetInput(control_type[i]); _assert_(input);
		}
		if(input->ObjectEnum()!=ControlInputEnum) _error_("input " << EnumToStringx(control_type[i]) << " is not a ControlInput");

		((ControlInput*)input)->UpdateValue(scalar);
		((ControlInput*)input)->Constrain();
		if (save_parameter) ((ControlInput*)input)->SaveValue();

		if(control_type[i]==MaterialsRheologyBbarEnum){
			this->InputExtrude(MaterialsRheologyBEnum,-1);
		}
		else if(control_type[i]==DamageDbarEnum){
			this->InputExtrude(DamageDEnum,-1);
		}
	}

	/*Clean up and return*/
cleanup_and_return:
	xDelete<int>(control_type);
}
/*}}}*/
void       Penta::InputDepthAverageAtBase(int enum_type,int average_enum_type){/*{{{*/

	int  step,i;
	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  Helem_list[NUMVERTICES];
	IssmDouble  zeros_list[NUMVERTICES]={0.0};
	IssmDouble  p0top1_list[NUMVERTICES];
	Penta* penta=NULL;
	Input* original_input=NULL;
	Input* element_integrated_input=NULL;
	Input* total_integrated_input=NULL;
	Input* element_thickness_input=NULL;
	Input* total_thickness_input=NULL;
	Input* depth_averaged_input=NULL;

	/*recover parameters: */

	/*Are we on the base? If not, return*/
	if(!IsOnBase()) return;

	/*OK, we are on bed. Initialize global inputs as 0*/
	total_thickness_input =new PentaInput(ThicknessEnum,zeros_list,P1Enum);

	/*Now follow all the upper element from the base to the surface to integrate the input*/
	penta=this;
	step =0;
	for(;;){

		/*Step1: Get original input (to be depth avegaged): */
		original_input=(Input*)penta->inputs->GetInput(enum_type);
		if(!original_input) _error_("could not find input with enum " << EnumToStringx(enum_type));

		/*If first time, initialize total_integrated_input*/
		if (step==0){
			if (original_input->ObjectEnum()==PentaInputEnum)
			 total_integrated_input=new PentaInput(average_enum_type,zeros_list,P1Enum);
			else if (original_input->ObjectEnum()==ControlInputEnum)
			 total_integrated_input=new PentaInput(average_enum_type,zeros_list,P1Enum);
			else if (original_input->ObjectEnum()==DoubleInputEnum)
			 total_integrated_input=new DoubleInput(average_enum_type,0.0);
			else{
			 _error_("object " << EnumToStringx(original_input->ObjectEnum()) << " not supported yet");
			}
		}

		/*Step2: Create element thickness input*/
		::GetVerticesCoordinates(&xyz_list[0][0],penta->vertices,NUMVERTICES);
		for(i=0;i<3;i++){
			Helem_list[i]=xyz_list[i+3][2]-xyz_list[i][2];
			Helem_list[i+3]=Helem_list[i];
		}
		element_thickness_input=new PentaInput(ThicknessEnum,Helem_list,P1Enum);

		/*Step3: Vertically integrate A COPY of the original*/
		if(original_input->ObjectEnum()==PentaInputEnum){
			if(((PentaInput*)original_input)->interpolation_type==P0Enum){
				original_input->GetInputValue(&p0top1_list[i]);
				element_integrated_input= new  PentaInput(original_input->InstanceEnum(),p0top1_list,P1Enum);
				element_integrated_input->VerticallyIntegrate(element_thickness_input);
			}
			else{
				element_integrated_input= (Input*)original_input->copy();
				element_integrated_input->VerticallyIntegrate(element_thickness_input);
			}
		}
		else{
			element_integrated_input= (Input*)original_input->copy();
			element_integrated_input->VerticallyIntegrate(element_thickness_input);
		}

		/*Add contributions to global inputs*/
		total_integrated_input->AXPY(element_integrated_input,1.0);
		total_thickness_input ->AXPY(element_thickness_input,1.0);

		/*Clean up*/
		delete element_thickness_input;
		delete element_integrated_input;

		/*Stop if we have reached the surface, else, take upper penta*/
		if (penta->IsOnSurface()) break;

		/* get upper Penta*/
		penta=penta->GetUpperPenta();
		_assert_(penta->Id()!=this->id);

		/*increase couter*/
		step++;
	}

	/*OK, now we only need to divide the depth integrated input by the total thickness!*/
	depth_averaged_input=total_integrated_input->PointwiseDivide(total_thickness_input);
	depth_averaged_input->ChangeEnum(average_enum_type);

	/*Clean up*/
	delete total_thickness_input;
	delete total_integrated_input;

	/*Finally, add to inputs*/
	this->inputs->AddInput((Input*)depth_averaged_input);
}
/*}}}*/
void       Penta::InputExtrude(int enum_type,int start){/*{{{*/

	_assert_(start==-1 || start==+1);

	/*Are we on the the boundary we want to be?*/
	if(start==-1 && !IsOnBase())    return;
	if(start==+1 && !IsOnSurface()) return;

	/*Step1: Get and Extrude original input: */
	Input* base_input=(Input*)this->inputs->GetInput(enum_type);
	if(!base_input) _error_("could not find input with enum " << EnumToStringx(enum_type));
	base_input->Extrude(start);

	/*Stop if there is only one layer of element*/
	if(start==-1 && this->IsOnSurface()) return;
	if(start==+1 && this->IsOnBase())    return;

	/*Step 2: this input has been extruded for this element, now follow the upper element*/
	Penta* penta=this;
	for(;;){
		/*get upper/lower Penta*/
		if(start==-1) penta=penta->GetUpperPenta();
		else          penta=penta->GetLowerPenta();
		_assert_(penta->Id()!=this->id);

		/*Add input of the basal element to penta->inputs*/
		Input* copy=(Input*)base_input->copy();
		penta->inputs->AddInput((Input*)copy);

		/*Stop if we have reached the surface/base*/
		if(start==-1 && penta->IsOnSurface()) break;
		if(start==+1 && penta->IsOnBase())    break;
	}
}
/*}}}*/
void       Penta::InputScale(int enum_type,IssmDouble scale_factor){/*{{{*/

	Input* input=NULL;

	/*Make a copy of the original input: */
	input=(Input*)this->inputs->GetInput(enum_type);
	if(!input)_error_("could not find old input with enum: " << EnumToStringx(enum_type));

	/*Scale: */
	input->Scale(scale_factor);
}
/*}}}*/
void       Penta::InputUpdateFromIoModel(int index,IoModel* iomodel){ /*{{{*/

	/*Intermediaries*/
	int         i,j;
	int         penta_vertex_ids[NUMVERTICES];
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
		penta_vertex_ids[i]=iomodel->elements[NUMVERTICES*index+i]; //ids for vertices are in the elements array from Matlab
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
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.balancethickness.thickening_rate")[penta_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i]/yts;
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i]/yts;
						this->inputs->AddInput(new ControlInput(BalancethicknessThickeningRateEnum,PentaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case VxEnum:
					if (iomodel->Data("md.initialization.vx")){
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.initialization.vx")[penta_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i]/yts;
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i]/yts;
						this->inputs->AddInput(new ControlInput(VxEnum,PentaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case VyEnum:
					if (iomodel->Data("md.initialization.vy")){
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.initialization.vy")[penta_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i]/yts;
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i]/yts;
						this->inputs->AddInput(new ControlInput(VyEnum,PentaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case FrictionCoefficientEnum:
					if (iomodel->Data("md.friction.coefficient")){
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.friction.coefficient")[penta_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i];
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i];
						this->inputs->AddInput(new ControlInput(FrictionCoefficientEnum,PentaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case MaterialsRheologyBbarEnum:
					if(iomodel->Data("md.materials.rheology_B")){
						for(j=0;j<NUMVERTICES;j++) nodeinputs[j]=iomodel->Data("md.materials.rheology_B")[penta_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i];
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i];
						this->inputs->AddInput(new ControlInput(MaterialsRheologyBEnum,PentaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case DamageDbarEnum:
					if(iomodel->Data("md.damage.D")){
						for(j=0;j<NUMVERTICES;j++) nodeinputs[j]=iomodel->Data("md.damage.D")[penta_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i];
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(penta_vertex_ids[j]-1)*num_control_type+i];
						this->inputs->AddInput(new ControlInput(DamageDEnum,PentaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
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
	if(control_analysis && iomodel->Data("md.inversion.cost_functions_coefficients")) {

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
			for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.inversion.cost_functions_coefficients")[(penta_vertex_ids[j]-1)*num_responses+i];
			datasetinput->AddInput(new PentaInput(InversionCostFunctionsCoefficientsEnum,nodeinputs,P1Enum),cost_functions_enums[i]);
		}

		/*Add datasetinput to element inputs*/
		this->inputs->AddInput(datasetinput);

		/*Free resources*/
		for(int j=0;j<num_cost_functions;j++) xDelete<char>(cost_functions[j]);
		xDelete<char*>(cost_functions);
		xDelete<int>(cost_functions_enums);
	}
}
/*}}}*/
void       Penta::InputUpdateFromSolutionOneDof(IssmDouble* solution,int enum_type){/*{{{*/

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
	this->inputs->AddInput(new PentaInput(enum_type,values,P1Enum));

	/*Free ressources:*/
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);
}
/*}}}*/
void       Penta::InputUpdateFromSolutionOneDofCollapsed(IssmDouble* solution,int enum_type){/*{{{*/

	const int  numdof   = NDOF1*NUMVERTICES;
	const int  numdof2d = NDOF1*NUMVERTICES2D;

	IssmDouble  values[numdof];
	int*    doflist = NULL;
	Penta  *penta   = NULL;

	/*If not on bed, return*/
	if (!IsOnBase()) return;

	/*Get dof list: */
	GetDofList(&doflist,NoneApproximationEnum,GsetEnum);

	/*Use the dof list to index into the solution vector and extrude it */
	for(int i=0;i<numdof2d;i++){
		values[i]         =solution[doflist[i]];
		values[i+numdof2d]=values[i];
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*Start looping over all elements above current element and update all inputs*/
	penta=this;
	for(;;){
		/*Add input to the element: */
		penta->inputs->AddInput(new PentaInput(enum_type,values,P1Enum));

		/*Stop if we have reached the surface*/
		if (penta->IsOnSurface()) break;

		/* get upper Penta*/
		penta=penta->GetUpperPenta(); _assert_(penta->Id()!=this->id);
	}

	/*Free ressources:*/
	xDelete<int>(doflist);
}
/*}}}*/
void       Penta::InputUpdateFromVector(IssmDouble* vector, int name, int type){/*{{{*/

	const int   numdof         = NDOF1 *NUMVERTICES;
	int        *doflist        = NULL;
	IssmDouble  values[numdof];

	/*Check that name is an element input*/
	if(!IsInput(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	switch(type){
		case VertexPIdEnum: 
			for (int i=0;i<NUMVERTICES;i++){
				values[i]=vector[this->vertices[i]->Pid()];
			}
			/*update input*/
			this->inputs->AddInput(new PentaInput(name,values,P1Enum));
			return;

		case VertexSIdEnum: 
			for (int i=0;i<NUMVERTICES;i++){
				values[i]=vector[this->vertices[i]->Sid()];
			}
			/*update input*/
			this->inputs->AddInput(new PentaInput(name,values,P1Enum));
			return;

		case NodesEnum:
			/*Get dof list: */
			GetDofList(&doflist,NoneApproximationEnum,GsetEnum);

			/*Use the dof list to index into the vector: */
			for(int i=0;i<numdof;i++){
				values[i]=vector[doflist[i]];
				if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
				if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
			}
			/*Add input to the element: */
			this->inputs->AddInput(new PentaInput(name,values,P1Enum));

			/*Free ressources:*/
			xDelete<int>(doflist);
			return;

	  case NodeSIdEnum:
			for(int i=0;i<NUMVERTICES;i++){
				values[i]=vector[nodes[i]->Sid()];
				if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
				if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
			}
			/*Add input to the element: */
			this->inputs->AddInput(new PentaInput(name,values,P1Enum));

			/*Free ressources:*/
			xDelete<int>(doflist);
			return;

	  default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}
}
/*}}}*/
bool       Penta::IsIcefront(void){/*{{{*/

	bool isicefront;
	int i,nrice;
   IssmDouble ls[NUMVERTICES];

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&ls[0],MaskIceLevelsetEnum);

	/* If only one vertex has ice, there is an ice front here */
	isicefront=false;
	if(IsIceInElement()){
		nrice=0;       
		for(i=0;i<NUMVERTICES2D;i++)
			if(ls[i]<0.) nrice++;
		if(nrice==1) isicefront= true;
	}
	return isicefront;
}/*}}}*/
bool       Penta::IsNodeOnShelfFromFlags(IssmDouble* flags){/*{{{*/

	int  i;
	bool shelf=false;

	for(i=0;i<NUMVERTICES;i++){
		if (flags[vertices[i]->Pid()]<0.){
			shelf=true;
			break;
		}
	}
	return shelf;
}
/*}}}*/
bool       Penta::IsOnBase(void){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);
	sum = values[0]+values[1]+values[2]+values[3]+values[4]+values[5];
	_assert_(sum==0. || sum==3.);

	if(sum==3){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
bool       Penta::IsOnSurface(void){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);
	sum = values[0]+values[1]+values[2]+values[3]+values[4]+values[5];
	_assert_(sum==0. || sum==3.);

	if(sum==3){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
bool       Penta::IsZeroLevelset(int levelset_enum){/*{{{*/

	bool        iszerols;
	IssmDouble  ls[NUMVERTICES];

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&ls[0],levelset_enum);

	/*If the level set has always same sign, there is no ice front here*/
	iszerols = false;
	if(IsIceInElement()){
		if(ls[0]*ls[1]<0. || ls[0]*ls[2]<0. || (ls[0]*ls[1]*ls[2]==0. && ls[0]*ls[1]+ls[0]*ls[2]+ls[1]*ls[2]<=0.)){
			iszerols = true;
		}
	}
	return iszerols;
}
/*}}}*/
void       Penta::JacobianDeterminant(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetJacobianDeterminant(pJdet,xyz_list,(GaussPenta*)gauss);

}
/*}}}*/
void       Penta::JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetTriaJacobianDeterminant(pJdet,xyz_list_base,(GaussPenta*)gauss);

}
/*}}}*/
void       Penta::JacobianDeterminantLine(IssmDouble* pJdet,IssmDouble* xyz_list_line,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetSegmentJacobianDeterminant(pJdet,xyz_list_line,(GaussPenta*)gauss);

}
/*}}}*/
void       Penta::JacobianDeterminantSurface(IssmDouble* pJdet,IssmDouble* xyz_list_quad,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetQuadJacobianDeterminant(pJdet,xyz_list_quad,(GaussPenta*)gauss);

}
/*}}}*/
void       Penta::JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_top,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetTriaJacobianDeterminant(pJdet,xyz_list_top,(GaussPenta*)gauss);

}
/*}}}*/
IssmDouble Penta::MassFlux(IssmDouble* segment){/*{{{*/

	IssmDouble mass_flux=0;

	if(!IsOnBase()) return mass_flux;

	/*Depth Averaging Vx and Vy*/
	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Spawn Tria element from the base of the Penta: */
	Tria* tria=(Tria*)SpawnTria(0,1,2);
	mass_flux=tria->MassFlux(segment);
	delete tria->material; delete tria;

	/*Delete Vx and Vy averaged*/
	this->inputs->DeleteInput(VxAverageEnum);
	this->inputs->DeleteInput(VyAverageEnum);

	/*clean up and return*/
	return mass_flux;
}
/*}}}*/
IssmDouble Penta::MassFlux(IssmDouble x1, IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id){/*{{{*/

	IssmDouble mass_flux=0;

	if(!IsOnBase()) return mass_flux;

	/*Depth Averaging Vx and Vy*/
	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Spawn Tria element from the base of the Penta: */
	Tria* tria=(Tria*)SpawnTria(0,1,2);
	mass_flux=tria->MassFlux(x1,y1,x2,y2,segment_id);
	delete tria->material; delete tria;

	/*Delete Vx and Vy averaged*/
	this->inputs->DeleteInput(VxAverageEnum);
	this->inputs->DeleteInput(VyAverageEnum);

	/*clean up and return*/
	return mass_flux;
}
/*}}}*/
IssmDouble Penta::MinEdgeLength(IssmDouble* xyz_list){/*{{{*/
	/*Return the minimum lenght of the nine egdes of the penta*/

	int    i,node0,node1;
	int    edges[9][2]={{0,1},{0,2},{1,2},{3,4},{3,5},{4,5},{0,3},{1,4},{2,5}}; //list of the nine edges
	IssmDouble length;
	IssmDouble minlength=-1;

	for(i=0;i<9;i++){
		/*Find the two nodes for this edge*/
		node0=edges[i][0];
		node1=edges[i][1];

		/*Compute the length of this edge and compare it to the minimal length*/
		length=sqrt(pow(xyz_list[node0*3+0]-xyz_list[node1*3+0],2)+pow(xyz_list[node0*3+1]-xyz_list[node1*3+1],2)+pow(xyz_list[node0*3+2]-xyz_list[node1*3+2],2));
		if(length<minlength || minlength<0) minlength=length;
	}

	return minlength;
}
/*}}}*/
Gauss*     Penta::NewGauss(void){/*{{{*/
	return new GaussPenta();
}
/*}}}*/
Gauss*     Penta::NewGauss(int order){/*{{{*/
	return new GaussPenta(order,order);
}
/*}}}*/
Gauss*     Penta::NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert){/*{{{*/

	IssmDouble  area_coordinates[4][3];

	GetAreaCoordinates(&area_coordinates[0][0],xyz_list_front,xyz_list,4);

	return new GaussPenta(area_coordinates,order_horiz,order_vert);
}
/*}}}*/
Gauss*     Penta::NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order){/*{{{*/
	return new GaussPenta(point1,fraction1,fraction2,mainlyfloating,order);
}
/*}}}*/
Gauss*     Penta::NewGaussBase(int order){/*{{{*/
	return new GaussPenta(0,1,2,order);
}
/*}}}*/
Gauss*     Penta::NewGaussLine(int vertex1,int vertex2,int order){/*{{{*/
	return new GaussPenta(vertex1,vertex2,order);
}
/*}}}*/
Gauss*     Penta::NewGaussTop(int order){/*{{{*/
	return new GaussPenta(3,4,5,order);
}
/*}}}*/
void       Penta::NodalFunctions(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,this->element_type);

}
/*}}}*/
void       Penta::NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussPenta*)gauss,this->element_type);

}
/*}}}*/
void       Penta::NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussPenta*)gauss,this->VelocityInterpolation());

}
/*}}}*/
void       Penta::NodalFunctionsMINIDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussPenta*)gauss,P1bubbleEnum);

}
/*}}}*/
void       Penta::NodalFunctionsPressure(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,this->PressureInterpolation());

}
/*}}}*/
void       Penta::NodalFunctionsP1(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,P1Enum);

}
/*}}}*/
void       Penta::NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussPenta*)gauss,P1Enum);

}
/*}}}*/
void       Penta::NodalFunctionsP2(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,P2Enum);

}
/*}}}*/
void       Penta::NodalFunctionsVelocity(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,this->VelocityInterpolation());

}
/*}}}*/
void       Penta::NodalFunctionsTensor(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,this->TensorInterpolation());

}
/*}}}*/
int        Penta::NodalValue(IssmDouble* pvalue, int index, int natureofdataenum){/*{{{*/

	int i;
	int found=0;
	IssmDouble value;
	Input* data=NULL;
	GaussPenta* gauss=NULL;

	/*First, serarch the input: */
	data=inputs->GetInput(natureofdataenum); 

	/*figure out if we have the vertex id: */
	found=0;
	for(i=0;i<NUMVERTICES;i++){
		if(index==vertices[i]->Id()){
			/*Do we have natureofdataenum in our inputs? :*/
			if(data){
				/*ok, we are good. retrieve value of input at vertex :*/
				gauss=new GaussPenta(); gauss->GaussVertex(i);
				data->GetInputValue(&value,gauss);
				found=1;
				break;
			}
		}
	}

	delete gauss;
	if(found)*pvalue=value;
	return found;
}
/*}}}*/
void       Penta::NormalBase(IssmDouble* bed_normal,IssmDouble* xyz_list){/*{{{*/

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
}
/*}}}*/
void       Penta::NormalSection(IssmDouble* normal,IssmDouble* xyz_list){/*{{{*/

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
void       Penta::NormalTop(IssmDouble* top_normal,IssmDouble* xyz_list){/*{{{*/

	int i;
	IssmDouble v13[3],v23[3];
	IssmDouble normal[3];
	IssmDouble normal_norm;

	for (i=0;i<3;i++){
		v13[i]=xyz_list[0*3+i]-xyz_list[2*3+i];
		v23[i]=xyz_list[1*3+i]-xyz_list[2*3+i];
	}

	normal[0]=v13[1]*v23[2]-v13[2]*v23[1];
	normal[1]=v13[2]*v23[0]-v13[0]*v23[2];
	normal[2]=v13[0]*v23[1]-v13[1]*v23[0];
	normal_norm=sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

	top_normal[0]=normal[0]/normal_norm;
	top_normal[1]=normal[1]/normal_norm;
	top_normal[2]=normal[2]/normal_norm;
}
/*}}}*/
int        Penta::NumberofNodesPressure(void){/*{{{*/
	return PentaRef::NumberofNodes(this->PressureInterpolation());
}
/*}}}*/
int        Penta::NumberofNodesVelocity(void){/*{{{*/
	return PentaRef::NumberofNodes(this->VelocityInterpolation());
}
/*}}}*/
int        Penta::ObjectEnum(void){/*{{{*/

	return PentaEnum;

}
/*}}}*/
void       Penta::PotentialUngrounding(Vector<IssmDouble>* potential_ungrounding){/*{{{*/

	IssmDouble  h[NUMVERTICES],r[NUMVERTICES],gl[NUMVERTICES];
	IssmDouble  bed_hydro;
	IssmDouble  rho_water,rho_ice,density;

	/*material parameters: */
	rho_water=matpar->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	density=rho_ice/rho_water;
	GetInputListOnVertices(&h[0],ThicknessEnum);
	GetInputListOnVertices(&r[0],BedEnum);
	GetInputListOnVertices(&gl[0],MaskGroundediceLevelsetEnum);

	/*go through vertices, and figure out which ones are on the ice sheet, and want to unground: */
	for(int i=0;i<NUMVERTICES;i++){
		/*Find if grounded vertices want to start floating*/
		if (gl[i]>0.){
			bed_hydro=-density*h[i];
			if(bed_hydro>r[i]){
				/*Vertex that could potentially unground, flag it*/
				potential_ungrounding->SetValue(vertices[i]->Pid(),1,INS_VAL);
			}
		}
	}
}
/*}}}*/
int        Penta::PressureInterpolation(void){/*{{{*/
	return PentaRef::PressureInterpolation(this->element_type);
}
/*}}}*/
void       Penta::ReduceMatrices(ElementMatrix* Ke,ElementVector* pe){/*{{{*/

	int analysis_type;
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	if(pe){
		if(analysis_type==StressbalanceAnalysisEnum){
			if(this->element_type==MINIcondensedEnum){
				int approximation;
				inputs->GetInputValue(&approximation,ApproximationEnum);
				if(approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
					//Do nothing, condensation already done in PVectorCoupling
				}
				else{
					int indices[3]={18,19,20};
					pe->StaticCondensation(Ke,3,&indices[0]);
				}
			}
			else if(this->element_type==P1bubblecondensedEnum){
				int size   = nodes[6]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
				int offset = 0;
				for(int i=0;i<6;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
				int* indices=xNew<int>(size);
				for(int i=0;i<size;i++) indices[i] = offset+i;
				pe->StaticCondensation(Ke,size,indices);
				xDelete<int>(indices);
			}
		}
	}

	if(Ke){
		if(analysis_type==StressbalanceAnalysisEnum){
			int approximation;
			inputs->GetInputValue(&approximation,ApproximationEnum);
			if(approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
				//Do nothing condensatino already done for Stokes part
			}
			else{
				if(this->element_type==MINIcondensedEnum){
					int indices[3]={18,19,20};
					Ke->StaticCondensation(3,&indices[0]);
				}
				else if(this->element_type==P1bubblecondensedEnum){
					int size   = nodes[6]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
					int offset = 0;
					for(int i=0;i<6;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
					int* indices=xNew<int>(size);
					for(int i=0;i<size;i++) indices[i] = offset+i;
					Ke->StaticCondensation(size,indices);
					xDelete<int>(indices);
				}
			}
		}
	}
}
/*}}}*/
void       Penta::ResetFSBasalBoundaryCondition(void){/*{{{*/

	int          approximation;
	int          numindices;
	int         *indices = NULL;
	IssmDouble   slopex,slopey,groundedice;
	IssmDouble   xz_plane[6];
	IssmDouble*  vertexapproximation= NULL;

	/*For FS only: we want the CS to be tangential to the bedrock*/
	inputs->GetInputValue(&approximation,ApproximationEnum);
	if(!IsOnBase() || (approximation!=FSApproximationEnum && approximation!=SSAFSApproximationEnum &&  approximation!=HOFSApproximationEnum)) return;

	/*Get number of nodes for velocity only and base*/
	BasalNodeIndices(&numindices,&indices,this->VelocityInterpolation());

	/*Get inputs*/
	Input* slopex_input=inputs->GetInput(BedSlopeXEnum); _assert_(slopex_input);
	Input* slopey_input=inputs->GetInput(BedSlopeYEnum); _assert_(slopey_input);
	Input* groundedicelevelset_input=inputs->GetInput(MaskGroundediceLevelsetEnum); _assert_(groundedicelevelset_input);

	/*Loop over basal nodes and update their CS*/
	GaussPenta* gauss = new GaussPenta();
	for(int i=0;i<numindices;i++){//FIXME
		gauss->GaussNode(this->VelocityInterpolation(),indices[i]);

		slopex_input->GetInputValue(&slopex,gauss);
		slopey_input->GetInputValue(&slopey,gauss);
		groundedicelevelset_input->GetInputValue(&groundedice,gauss);

		/*New X axis          New Z axis*/
		xz_plane[0]=1.;       xz_plane[3]=-slopex;  
		xz_plane[1]=0.;       xz_plane[4]=-slopey;  
		xz_plane[2]=slopex;   xz_plane[5]=1.;          

		if(groundedice>=0){
			if(this->nodes[indices[i]]->GetApproximation()==FSvelocityEnum){
				this->nodes[indices[i]]->DofInSSet(2); //vz 
			}
			else if(this->nodes[indices[i]]->GetApproximation()==SSAFSApproximationEnum || this->nodes[indices[i]]->GetApproximation()==HOFSApproximationEnum){
				this->nodes[indices[i]]->DofInSSet(4); //vz 
			}
			else _error_("Flow equation approximation"<<EnumToStringx(this->nodes[indices[i]]->GetApproximation())<<" not supported yet");
		}
		else{
			if(this->nodes[indices[i]]->GetApproximation()==FSvelocityEnum){
				this->nodes[indices[i]]->DofInFSet(2); //vz
			}
			else if(this->nodes[indices[i]]->GetApproximation()==SSAFSApproximationEnum || this->nodes[indices[i]]->GetApproximation()==HOFSApproximationEnum){
				this->nodes[indices[i]]->DofInFSet(4); //vz 
			}
			else _error_("Flow equation approximation"<<EnumToStringx(this->nodes[indices[i]]->GetApproximation())<<" not supported yet");
		}

		XZvectorsToCoordinateSystem(&this->nodes[indices[i]]->coord_system[0][0],&xz_plane[0]);
	}

	/*cleanup*/
	xDelete<int>(indices);
	delete gauss;
}
/*}}}*/
void       Penta::ResetHooks(){/*{{{*/

	this->nodes=NULL;
	this->vertices=NULL;
	this->material=NULL;
	this->matpar=NULL;
	this->verticalneighbors=NULL;
	this->parameters=NULL;

	//deal with ElementHook mother class
	for(int i=0;i<this->numanalyses;i++) if(this->hnodes[i]) this->hnodes[i]->reset();
	this->hvertices->reset();
	this->hmaterial->reset();
	this->hmatpar->reset();
	if(this->hneighbors) this->hneighbors->reset();

}
/*}}}*/
void       Penta::ResetLevelsetFromSegmentlist(IssmDouble* segments,int numsegments){/*{{{*/

	/*Intermediaries*/
	IssmDouble d,xn,yn;

	/*Get current levelset and vertex coordinates*/
	IssmDouble ls[NUMVERTICES];
	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	GetInputListOnVertices(&ls[0],MaskIceLevelsetEnum);
	InputDuplicate(MaskIceLevelsetEnum,PressureEnum);

	/*Get distance from list of segments and reset ls*/
	for(int j=0;j<NUMVERTICES;j++){
		IssmDouble dmin = 1.e+50;
		for(int i=0;i<numsegments;i++){
			IssmDouble x = xyz_list[j][0];
			IssmDouble y = xyz_list[j][1];
			IssmDouble l2 = (segments[4*i+2]-segments[4*i+0])*(segments[4*i+2]-segments[4*i+0]) + (segments[4*i+3]-segments[4*i+1])*(segments[4*i+3]-segments[4*i+1]);

			/*Segment has a length of 0*/
			if(l2==0.){
				d = (x-segments[4*i+0])*(x-segments[4*i+0])+(y-segments[4*i+1])*(y-segments[4*i+1]);
				if(d<dmin) dmin = d;
				continue;
			}

			/*Consider the line extending the segment, parameterized as v + t (w - v).
			 *We find projection of point p onto the line.
			 *It falls where t = [(p-v) . (w-v)] / |w-v|^2*/
			IssmDouble t = ((x-segments[4*i+0])*(segments[4*i+2]-segments[4*i+0]) + (y-segments[4*i+1])*(segments[4*i+3]-segments[4*i+1]))/l2;
			if(t < 0.0){
				// Beyond the 'v' end of the segment
				d = (x-segments[4*i+0])*(x-segments[4*i+0])+(y-segments[4*i+1])*(y-segments[4*i+1]);
			}
			else if (t > 1.0){
				// Beyond the 'w' end of the segment
				d = (x-segments[4*i+2])*(x-segments[4*i+2])+(y-segments[4*i+3])*(y-segments[4*i+3]);
			}
			else{
				// Projection falls on the segment
				xn = segments[4*i+0] + t * (segments[4*i+2] - segments[4*i+0]);
				yn = segments[4*i+1] + t * (segments[4*i+3] - segments[4*i+1]);
				d = (x-xn)*(x-xn)+(y-yn)*(y-yn);
			}

			if(d<dmin) dmin = d;
		}

		/*Update signed distance*/
		dmin = sqrt(dmin);
		if(dmin>10000) dmin=10000;
		if(ls[j]>0){
			ls[j] = dmin;
		}
		else{
			ls[j] = - dmin;
		}
	}

	/*Update Levelset*/
	this->inputs->AddInput(new PentaInput(MaskIceLevelsetEnum,&ls[0],P1Enum));
}
/*}}}*/
void       Penta::SetClone(int* minranks){/*{{{*/

	_error_("not implemented yet");
}

/*}}}*/
void       Penta::SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index){/*{{{*/

	IssmDouble  values[NUMVERTICES];
	int         vertexpidlist[NUMVERTICES],control_init;

	/*Specific case for depth averaged quantities*/
	control_init=control_enum;
	if(control_enum==MaterialsRheologyBbarEnum){
		control_enum=MaterialsRheologyBEnum;
		if(!IsOnBase()) return;
	}
	if(control_enum==DamageDbarEnum){
		control_enum=DamageDEnum;
		if(!IsOnBase()) return;
	}

	/*Get out if this is not an element input*/
	if(!IsInput(control_enum)) return;

	/*Prepare index list*/
	GradientIndexing(&vertexpidlist[0],control_index);

	/*Get values on vertices*/
	for(int i=0;i<NUMVERTICES;i++){
		values[i]=vector[vertexpidlist[i]];
	}
	Input* new_input = new PentaInput(control_enum,values,P1Enum);
	Input* input=(Input*)this->inputs->GetInput(control_enum);   _assert_(input);
	if(input->ObjectEnum()!=ControlInputEnum){
		_error_("input " << EnumToStringx(control_enum) << " is not a ControlInput");
	}

	((ControlInput*)input)->SetInput(new_input);

	if(control_init==MaterialsRheologyBbarEnum){
		this->InputExtrude(control_enum,-1);
	}
	if(control_init==DamageDbarEnum){
		this->InputExtrude(control_enum,-1);
	}
}
/*}}}*/
void       Penta::SetCurrentConfiguration(Elements* elementsin, Loads* loadsin, Nodes* nodesin, Materials* materialsin, Parameters* parametersin){/*{{{*/

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
void       Penta::SetTemporaryElementType(int element_type_in){/*{{{*/
	this->element_type=element_type_in;
}
/*}}}*/
Element*   Penta::SpawnBasalElement(void){/*{{{*/

	_assert_(this->IsOnBase());

	switch(this->material->ObjectEnum()){
		case MaticeEnum:
			this->InputDepthAverageAtBase(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
			if(this->material->IsDamage())this->InputDepthAverageAtBase(DamageDEnum,DamageDbarEnum);
			if(this->material->IsEnhanced())this->InputDepthAverageAtBase(MaterialsRheologyEEnum,MaterialsRheologyEbarEnum);
			break;
		case MatestarEnum:
			this->InputDepthAverageAtBase(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
			this->InputDepthAverageAtBase(MaterialsRheologyEcEnum,MaterialsRheologyEcbarEnum);
			this->InputDepthAverageAtBase(MaterialsRheologyEsEnum,MaterialsRheologyEsbarEnum);
			break;
		default:
			_error_("not supported yet");
	}
	if(this->inputs->GetInput(VxEnum)) this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	if(this->inputs->GetInput(VyEnum)) this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);
	if(this->inputs->GetInput(CalvingratexEnum)) this->InputDepthAverageAtBase(CalvingratexEnum,CalvingratexAverageEnum);
	if(this->inputs->GetInput(CalvingrateyEnum)) this->InputDepthAverageAtBase(CalvingrateyEnum,CalvingrateyAverageEnum);
	Tria* tria=(Tria*)SpawnTria(0,1,2);
	switch(this->material->ObjectEnum()){
		case MaticeEnum:
			this->inputs->DeleteInput(MaterialsRheologyBbarEnum);
			this->inputs->DeleteInput(DamageDbarEnum);
			break;
		case MatestarEnum:
			break;
		default:
			_error_("not supported yet");
	}
	this->inputs->DeleteInput(VxAverageEnum);
	this->inputs->DeleteInput(VyAverageEnum);
	this->inputs->DeleteInput(CalvingratexAverageEnum);
	this->inputs->DeleteInput(CalvingrateyAverageEnum);

	return tria;
}
/*}}}*/
Element*   Penta::SpawnTopElement(void){/*{{{*/

	_assert_(this->IsOnSurface());

	Tria* tria=(Tria*)SpawnTria(3,4,5);

	return tria;
}
/*}}}*/
Tria*      Penta::SpawnTria(int index1,int index2,int index3){/*{{{*/

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
	tria->nodes=(Node**)tria->hnodes[analysis_counter]->deliverp();
	tria->vertices=(Vertex**)tria->hvertices->deliverp();
	tria->matpar=(Matpar*)tria->hmatpar->delivers();

	/*Return new Tria*/
	return tria;
}
/*}}}*/
IssmDouble Penta::StabilizationParameter(IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble diameter, IssmDouble kappa){/*{{{*/
	/*Compute stabilization parameter*/
	/*kappa=thermalconductivity/(rho_ice*hearcapacity) for thermal model*/
	/*kappa=enthalpydiffusionparameter for enthalpy model*/

	IssmDouble normu;
	IssmDouble tau_parameter;

	normu=pow(pow(u,2)+pow(v,2)+pow(w,2),0.5);
	if(normu*diameter/(3*2*kappa)<1){ 
		tau_parameter=pow(diameter,2)/(3*2*2*kappa);
	}
	else tau_parameter=diameter/(2*normu);

	return tau_parameter;
}
/*}}}*/
void       Penta::StrainRateparallel(){/*{{{*/

	IssmDouble *xyz_list = NULL;
	IssmDouble  epsilon[6];
	GaussPenta* gauss=NULL;
	IssmDouble  vx,vy,vel;
	IssmDouble  strainxx;
	IssmDouble  strainxy;
	IssmDouble  strainyy;
	IssmDouble  strainparallel[NUMVERTICES];

	/* Get node coordinates and dof list: */
	this->GetVerticesCoordinates(&xyz_list);

	/*Retrieve all inputs we will need*/
	Input* vx_input=inputs->GetInput(VxEnum);                                  _assert_(vx_input);
	Input* vy_input=inputs->GetInput(VyEnum);                                  _assert_(vy_input);
	Input* vz_input=inputs->GetInput(VzEnum);												_assert_(vz_input);

	/* Start looping on the number of vertices: */
	gauss=new GaussPenta();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/* Get the value we need*/
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vel=vx*vx+vy*vy;

		/*Compute strain rate and viscosity: */
		this->StrainRateFS(&epsilon[0],xyz_list,gauss,vx_input,vy_input,vz_input);
		strainxx=epsilon[0];
		strainyy=epsilon[1];
		strainxy=epsilon[3];

		/*strainparallel= Strain rate along the ice flow direction */
		strainparallel[iv]=(vx*vx*(strainxx)+vy*vy*(strainyy)+2*vy*vx*strainxy)/(vel+1.e-14);
	}

	/*Add input*/
	this->inputs->AddInput(new PentaInput(StrainRateparallelEnum,&strainparallel[0],P1Enum));

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
}
/*}}}*/
void       Penta::StrainRateperpendicular(){/*{{{*/

	IssmDouble *xyz_list = NULL;
	IssmDouble  epsilon[6];
	GaussPenta* gauss=NULL;
	IssmDouble  vx,vy,vel;
	IssmDouble  strainxx;
	IssmDouble  strainxy;
	IssmDouble  strainyy;
	IssmDouble  strainperpendicular[NUMVERTICES];

	/* Get node coordinates and dof list: */
	this->GetVerticesCoordinates(&xyz_list);

	/*Retrieve all inputs we will need*/
	Input* vx_input=inputs->GetInput(VxEnum);                                  _assert_(vx_input);
	Input* vy_input=inputs->GetInput(VyEnum);                                  _assert_(vy_input);
	Input* vz_input=inputs->GetInput(VzEnum);												_assert_(vz_input);

	/* Start looping on the number of vertices: */
	gauss=new GaussPenta();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/* Get the value we need*/
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vel=vx*vx+vy*vy;

		/*Compute strain rate and viscosity: */
		this->StrainRateFS(&epsilon[0],xyz_list,gauss,vx_input,vy_input,vz_input);
		strainxx=epsilon[0];
		strainyy=epsilon[1];
		strainxy=epsilon[3];

		/*strainperpendicular= Strain rate perpendicular to the ice flow direction */
		strainperpendicular[iv]=(vx*vx*(strainyy)+vy*vy*(strainxx)-2*vy*vx*strainxy)/(vel+1.e-14);
	}

	/*Add input*/
	this->inputs->AddInput(new PentaInput(StrainRateperpendicularEnum,&strainperpendicular[0],P1Enum));

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
}
/*}}}*/
void       Penta::StressIntensityFactor(){/*{{{*/

	/* Check if we are on the base */
	if(!IsOnBase()) return;

	IssmDouble  ki[6]={0.};
	IssmDouble  const_grav=9.81;
	IssmDouble  rho_ice=900;
	IssmDouble  rho_water=1000;
	IssmDouble  Jdet[3];
	IssmDouble  pressure,vx,vy,vel,deviaxx,deviaxy,deviayy,water_depth,prof,stress_xx,thickness;

	Penta* penta=this;
	for(;;){
	
		IssmDouble  xyz_list[NUMVERTICES][3];
		/* Get node coordinates and dof list: */
		::GetVerticesCoordinates(&xyz_list[0][0],penta->vertices,NUMVERTICES);

		///*Compute the Jacobian for the vertical integration*/
		Jdet[0]=(xyz_list[3][2]-xyz_list[0][2])*0.5;
		Jdet[1]=(xyz_list[4][2]-xyz_list[1][2])*0.5;
		Jdet[2]=(xyz_list[5][2]-xyz_list[2][2])*0.5;
	
		/*Retrieve all inputs we will need*/
		Input* vx_input=inputs->GetInput(VxEnum);                                  _assert_(vx_input);
		Input* vy_input=inputs->GetInput(VyEnum);                                  _assert_(vy_input);
		Input* vel_input=inputs->GetInput(VelEnum);                                _assert_(vel_input);
		Input* pressure_input=inputs->GetInput(PressureEnum);                      _assert_(pressure_input);
		Input* deviaxx_input=inputs->GetInput(DeviatoricStressxxEnum);             _assert_(deviaxx_input);
		Input* deviaxy_input=inputs->GetInput(DeviatoricStressxyEnum);             _assert_(deviaxy_input);
		Input* deviayy_input=inputs->GetInput(DeviatoricStressyyEnum);             _assert_(deviayy_input);
		Input* surface_input=inputs->GetInput(SurfaceEnum);								_assert_(surface_input);
		Input* thickness_input=inputs->GetInput(ThicknessEnum);							_assert_(thickness_input);
		
		/* Start looping on the number of 2D vertices: */
		for(int ig=0;ig<3;ig++){
			GaussPenta* gauss=new GaussPenta(ig,3+ig,11);
			for (int iv=gauss->begin();iv<gauss->end();iv++){
				gauss->GaussPoint(iv);

				/* Get the value we need*/
				pressure_input->GetInputValue(&pressure,gauss);
				vx_input->GetInputValue(&vx,gauss);
				vy_input->GetInputValue(&vy,gauss);
				vel_input->GetInputValue(&vel,gauss);
				deviaxx_input->GetInputValue(&deviaxx,gauss);
				deviaxy_input->GetInputValue(&deviaxy,gauss);
				deviayy_input->GetInputValue(&deviayy,gauss);
				surface_input->GetInputValue(&water_depth,gauss);
				thickness_input->GetInputValue(&thickness,gauss);
				prof=water_depth-penta->GetZcoord(&xyz_list[0][0],gauss);

				/*stress_xx= Deviatoric stress along the ice flow direction plus cryostatic pressure */
				stress_xx=(vx*vx*(deviaxx)+vy*vy*(deviayy)+2*vy*vx*deviaxy)/(vel*vel+1.e-6);

				if(prof<water_depth&prof<thickness){
					/* Compute the local stress intensity factor*/ 
					ki[ig]+=Jdet[ig]*gauss->weight*stress_xx*StressIntensityIntegralWeight(prof,min(water_depth,thickness),thickness);
				}
			}
			delete gauss;
		}
			
		/*Stop if we have reached the surface/base*/
		if(penta->IsOnSurface()) break;
		
		/*get upper Penta*/
		penta=penta->GetUpperPenta();
		_assert_(penta->Id()!=this->id);
	}

	/*Add input*/
	this->inputs->AddInput(new PentaInput(StressIntensityFactorEnum,&ki[0],P1Enum));
	this->InputExtrude(StressIntensityFactorEnum,-1);
}
/*}}}*/
IssmDouble Penta::SurfaceArea(void){/*{{{*/

	int    approximation;
	IssmDouble S;
	Tria*  tria=NULL;

	/*retrieve inputs :*/
	inputs->GetInputValue(&approximation,ApproximationEnum);

	/*If on water, return 0: */
	if(!IsIceInElement())return 0;

	/*Bail out if this element if:
	 * -> Non SSA not on the surface
	 * -> SSA (2d model) and not on bed) */
	if ((approximation!=SSAApproximationEnum && !IsOnSurface()) || (approximation==SSAApproximationEnum && !IsOnBase())){
		return 0;
	}
	else if (approximation==SSAApproximationEnum){

		/*This element should be collapsed into a tria element at its base. Create this tria element, 
		 * and compute SurfaceArea*/
		tria=(Tria*)SpawnTria(0,1,2);
		S=tria->SurfaceArea();
		delete tria->material; delete tria;
		return S;
	}
	else{

		tria=(Tria*)SpawnTria(3,4,5);
		S=tria->SurfaceArea();
		delete tria->material; delete tria;
		return S;
	}
}
/*}}}*/
IssmDouble Penta::TimeAdapt(void){/*{{{*/

	int    i;
	IssmDouble C,dx,dy,dz,dt;
	IssmDouble maxabsvx,maxabsvy,maxabsvz;
	IssmDouble maxx,minx,maxy,miny,maxz,minz;
	IssmDouble xyz_list[NUMVERTICES][3];

	/*get CFL coefficient:*/
	this->parameters->FindParam(&C,TimesteppingCflCoefficientEnum);

	/*Get for Vx and Vy, the max of abs value: */
	maxabsvx = this->inputs->MaxAbs(VxEnum);
	maxabsvy = this->inputs->MaxAbs(VyEnum);
	maxabsvz = this->inputs->MaxAbs(VzEnum);

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	minx=xyz_list[0][0];
	maxx=xyz_list[0][0];
	miny=xyz_list[0][1];
	maxy=xyz_list[0][1];
	minz=xyz_list[0][2];
	maxz=xyz_list[0][2];

	for(i=1;i<NUMVERTICES;i++){
		if (xyz_list[i][0]<minx)minx=xyz_list[i][0];
		if (xyz_list[i][0]>maxx)maxx=xyz_list[i][0];
		if (xyz_list[i][1]<miny)miny=xyz_list[i][1];
		if (xyz_list[i][1]>maxy)maxy=xyz_list[i][1];
		if (xyz_list[i][2]<minz)minz=xyz_list[i][2];
		if (xyz_list[i][2]>maxz)maxz=xyz_list[i][2];
	}
	dx=maxx-minx;
	dy=maxy-miny;
	dz=maxz-minz;

	/*CFL criterion: */
	dt=C/(maxabsvx/dx+maxabsvy/dy+maxabsvz/dz);

	return dt;
}/*}}}*/
IssmDouble Penta::TotalFloatingBmb(void){/*{{{*/

	/*The fbmb[kg yr-1] of one element is area[m2] * melting_rate [kg m^-2 yr^-1]*/
	int        point1;
	bool       mainlyfloating;
	IssmDouble fbmb=0;
	IssmDouble rho_ice,fraction1,fraction2,floatingmelt,Jdet;
	IssmDouble Total_Fbmb=0;
	IssmDouble xyz_list[NUMVERTICES][3];
	Gauss*     gauss     = NULL;

   if(!IsIceInElement() || !IsOnBase())return 0;

	/*Get material parameters :*/
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	Input* floatingmelt_input = this->GetInput(BasalforcingsFloatingiceMeltingRateEnum); _assert_(floatingmelt_input); 
	Input* gllevelset_input = this->GetInput(MaskGroundediceLevelsetEnum); _assert_(gllevelset_input);
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	this->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating);
	/* Start  looping on the number of gaussian points: */
	gauss = this->NewGauss(point1,fraction1,fraction2,1-mainlyfloating,3);
	for(int ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);
		this->JacobianDeterminantBase(&Jdet,&xyz_list[0][0],gauss);
		floatingmelt_input->GetInputValue(&floatingmelt,gauss);
		fbmb+=floatingmelt*Jdet*gauss->weight;
	}

   Total_Fbmb=rho_ice*fbmb;	        // from volume to mass

	/*Return: */
	delete gauss;
	return Total_Fbmb;
}
/*}}}*/
IssmDouble Penta::TotalGroundedBmb(void){/*{{{*/

	/*The gbmb[kg yr-1] of one element is area[m2] * gounded melting rate [kg m^-2 yr^-1]*/
	int        point1;
	bool       mainlyfloating;
	IssmDouble gbmb=0;
	IssmDouble rho_ice,fraction1,fraction2,groundedmelt,Jdet;
	IssmDouble Total_Gbmb=0;
	IssmDouble xyz_list[NUMVERTICES][3];
	Gauss*     gauss     = NULL;

   if(!IsIceInElement() || !IsOnBase())return 0;

	/*Get material parameters :*/
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	Input* groundedmelt_input = this->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(groundedmelt_input); 
	Input* gllevelset_input = this->GetInput(MaskGroundediceLevelsetEnum); _assert_(gllevelset_input);
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	this->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating);
	/* Start  looping on the number of gaussian points: */
	gauss = this->NewGauss(point1,fraction1,fraction2,mainlyfloating,3);
	for(int ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);
		this->JacobianDeterminantBase(&Jdet,&xyz_list[0][0],gauss);
		groundedmelt_input->GetInputValue(&groundedmelt,gauss);
		gbmb+=groundedmelt*Jdet*gauss->weight;
	}

   Total_Gbmb=rho_ice*gbmb;	        // from volume to mass

	/*Return: */
	delete gauss;
	return Total_Gbmb;
}
/*}}}*/
IssmDouble Penta::TotalSmb(void){/*{{{*/

	/*The smb[Gt yr-1] of one element is area[m2] * smb [ m ice yr^-1] * rho_ice [kg m-3] / 1e+10^12 */
	IssmDouble base,smb,rho_ice;
	IssmDouble Total_Smb=0;
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Get material parameters :*/
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);

	if(!IsIceInElement() || !IsOnSurface()) return 0.;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Triangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2. * fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	/*Now get the average SMB over the element*/
	Input* smb_input = inputs->GetInput(SmbMassBalanceEnum); _assert_(smb_input);

	smb_input->GetInputAverage(&smb);
	Total_Smb=rho_ice*base*smb;// smb on element in kg s-1

	/*Return: */
	return Total_Smb;
}
/*}}}*/
void       Penta::Update(int index,IoModel* iomodel,int analysis_counter,int analysis_type,int finiteelement_type){ /*{{{*/

	/*Intermediaries*/
	int        i;
	int        penta_vertex_ids[6];
	IssmDouble nodeinputs[6];
	IssmDouble yts;
	bool       dakota_analysis;
	int        numnodes;
	int*       penta_node_ids = NULL;

	/*Fetch parameters: */
	iomodel->FindConstant(&yts,"md.constants.yts");
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");

	/*Checks if debuging*/
	_assert_(iomodel->elements);

	/*Recover element type*/
	this->element_type_list[analysis_counter]=finiteelement_type;

	/*Recover vertices ids needed to initialize inputs*/
	for(i=0;i<6;i++) penta_vertex_ids[i]=iomodel->elements[6*index+i]; //ids for vertices are in the elements array from Matlab

	/*Recover nodes ids needed to initialize the node hook.*/
	switch(finiteelement_type){
		case P1Enum:
			numnodes         = 6;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			break;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			numnodes         = 7;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[6]=iomodel->nodecounter+iomodel->numberofvertices+index+1;
			break;
		case P1xP2Enum:
			numnodes         = 9;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			break;
		case P1xP3Enum:
			numnodes         = 12;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+2*iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+2*iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+2*iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+2*iomodel->elementtoedgeconnectivity[9*index+0]+2;
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+2*iomodel->elementtoedgeconnectivity[9*index+1]+2;
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+2*iomodel->elementtoedgeconnectivity[9*index+2]+2;
			break;
		case P2xP1Enum:
			numnodes         = 12;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			break;
		case P2xP4Enum:
			numnodes         = 30;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0]; /*Vertex 1*/
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1]; /*Vertex 2*/
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2]; /*Vertex 3*/
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3]; /*Vertex 4*/
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4]; /*Vertex 5*/
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5]; /*Vertex 6*/
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+0]+1; /*mid vertical edge 1*/
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+1]+1; /*mid vertical edge 2*/
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+2]+1; /*mid vertical edge 3*/
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+3]+1; /*mid basal edge 1*/
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+4]+1; /*mid basal edge 2*/
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+5]+1; /*mid basal edge 3*/
			penta_node_ids[12]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+6]+1; /*mid top edge 1*/
			penta_node_ids[13]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+7]+1; /*mid top edge 2*/
			penta_node_ids[14]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+8]+1; /*mid top edge 3*/
			penta_node_ids[15]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+0]+2; /* 1/4 vertical edge 1*/
			penta_node_ids[16]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+1]+2; /* 1/4 vertical edge 2*/
			penta_node_ids[17]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+2]+2; /* 1/4 vertical edge 3*/
			penta_node_ids[18]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+0]+3; /* 2/4 vertical edge 1*/
			penta_node_ids[19]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+1]+3; /* 2/4 vertical edge 2*/
			penta_node_ids[20]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+2]+3; /* 2/4 vertical edge 3*/
			penta_node_ids[21]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+2]+1; /* 1/4 vertical face 1*/
			penta_node_ids[22]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+3]+1; /* 1/4 vertical face 2*/
			penta_node_ids[23]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+4]+1; /* 1/4 vertical face 3*/
			penta_node_ids[24]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+2]+2; /* 2/4 vertical face 1*/
			penta_node_ids[25]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+3]+2; /* 2/4 vertical face 2*/
			penta_node_ids[26]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+4]+2; /* 2/4 vertical face 3*/
			penta_node_ids[27]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+2]+3; /* 3/4 vertical face 1*/
			penta_node_ids[28]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+3]+3; /* 3/4 vertical face 2*/
			penta_node_ids[29]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+4]+3; /* 3/4 vertical face 3*/
			break;
		case P2Enum:
			numnodes         = 18;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+2]+1;
			penta_node_ids[16]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+3]+1;
			penta_node_ids[17]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+4]+1;
			break;
		case P2bubbleEnum: case P2bubblecondensedEnum:
			numnodes         = 19;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+2]+1;
			penta_node_ids[16]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+3]+1;
			penta_node_ids[17]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+4]+1;
			penta_node_ids[18]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+index+1;
			break;
		case P1P1Enum: case P1P1GLSEnum:
			numnodes         = 12;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];

			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elements[6*index+0];
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elements[6*index+1];
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elements[6*index+2];
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elements[6*index+3];
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elements[6*index+4];
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elements[6*index+5];
			break;
		case MINIEnum: case MINIcondensedEnum:
			numnodes         = 13;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+index+1;

			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+0];
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+1];
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+2];
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+3];
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+4];
			penta_node_ids[12]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+5];
			break;
		case TaylorHoodEnum:
			numnodes         = 24;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+2]+1;
			penta_node_ids[16]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+3]+1;
			penta_node_ids[17]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+4]+1;

			penta_node_ids[18]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->elements[6*index+0];
			penta_node_ids[19]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->elements[6*index+1];
			penta_node_ids[20]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->elements[6*index+2];
			penta_node_ids[21]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->elements[6*index+3];
			penta_node_ids[22]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->elements[6*index+4];
			penta_node_ids[23]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->elements[6*index+5];
			break;
		case LATaylorHoodEnum:
			numnodes         = 18;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+2]+1;
			penta_node_ids[16]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+3]+1;
			penta_node_ids[17]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+4]+1;
			break;
		case OneLayerP4zEnum:
			numnodes         = 30+6;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0]; /*Vertex 1*/
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1]; /*Vertex 2*/
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2]; /*Vertex 3*/
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3]; /*Vertex 4*/
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4]; /*Vertex 5*/
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5]; /*Vertex 6*/
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+0]+1; /*mid vertical edge 1*/
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+1]+1; /*mid vertical edge 2*/
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+2]+1; /*mid vertical edge 3*/
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+3]+1; /*mid basal edge 1*/
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+4]+1; /*mid basal edge 2*/
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+5]+1; /*mid basal edge 3*/
			penta_node_ids[12]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+6]+1; /*mid top edge 1*/
			penta_node_ids[13]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+7]+1; /*mid top edge 2*/
			penta_node_ids[14]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+8]+1; /*mid top edge 3*/
			penta_node_ids[15]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+0]+2; /* 1/4 vertical edge 1*/
			penta_node_ids[16]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+1]+2; /* 1/4 vertical edge 2*/
			penta_node_ids[17]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+2]+2; /* 1/4 vertical edge 3*/
			penta_node_ids[18]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+0]+3; /* 2/4 vertical edge 1*/
			penta_node_ids[19]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+1]+3; /* 2/4 vertical edge 2*/
			penta_node_ids[20]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->elementtoedgeconnectivity[9*index+2]+3; /* 2/4 vertical edge 3*/
			penta_node_ids[21]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+2]+1; /* 1/4 vertical face 1*/
			penta_node_ids[22]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+3]+1; /* 1/4 vertical face 2*/
			penta_node_ids[23]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+4]+1; /* 1/4 vertical face 3*/
			penta_node_ids[24]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+2]+2; /* 2/4 vertical face 1*/
			penta_node_ids[25]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+3]+2; /* 2/4 vertical face 2*/
			penta_node_ids[26]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+4]+2; /* 2/4 vertical face 3*/
			penta_node_ids[27]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+2]+3; /* 3/4 vertical face 1*/
			penta_node_ids[28]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+3]+3; /* 3/4 vertical face 2*/
			penta_node_ids[29]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->elementtofaceconnectivity[5*index+4]+3; /* 3/4 vertical face 3*/

			penta_node_ids[30]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->numberoffaces+iomodel->elements[6*index+0];
			penta_node_ids[31]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->numberoffaces+iomodel->elements[6*index+1];
			penta_node_ids[32]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->numberoffaces+iomodel->elements[6*index+2];
			penta_node_ids[33]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->numberoffaces+iomodel->elements[6*index+3];
			penta_node_ids[34]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->numberoffaces+iomodel->elements[6*index+4];
			penta_node_ids[35]=iomodel->nodecounter+iomodel->numberofvertices+3*iomodel->numberofedges+3*iomodel->numberoffaces+iomodel->elements[6*index+5];
			break;
		case CrouzeixRaviartEnum:
			numnodes         = 25;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+2]+1;
			penta_node_ids[16]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+3]+1;
			penta_node_ids[17]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+4]+1;
			penta_node_ids[18]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+index+1;

			penta_node_ids[19]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+1;
			penta_node_ids[20]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+2;
			penta_node_ids[21]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+3;
			penta_node_ids[22]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+4;
			penta_node_ids[23]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+5;
			penta_node_ids[24]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+6;
			break;
		case LACrouzeixRaviartEnum:
			numnodes         = 19;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->nodecounter+iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->nodecounter+iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->nodecounter+iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->nodecounter+iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->nodecounter+iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->nodecounter+iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+2]+1;
			penta_node_ids[16]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+3]+1;
			penta_node_ids[17]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+4]+1;
			penta_node_ids[18]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+index+1;
			break;
		default:
			_error_("Finite element "<<EnumToStringx(finiteelement_type)<<" not supported yet");
	}

	/*hooks: */
	this->SetHookNodes(penta_node_ids,numnodes,analysis_counter); this->nodes=NULL;
	xDelete<int>(penta_node_ids);

	/*Fill with IoModel*/
	this->InputUpdateFromIoModel(index,iomodel);

	/*Defaults if not provided in iomodel*/
	switch(analysis_type){

		case StressbalanceAnalysisEnum:
			_assert_(iomodel->Data("md.flowequation.element_equation"));

			if((IoCodeToEnumElementEquation(reCast<int>(iomodel->Data("md.flowequation.element_equation")[index])))==HOFSApproximationEnum){
				/*Create VzHO and VzFS Enums*/
				if(iomodel->Data("md.initialization.vz") && iomodel->Data("md.flowequation.borderFS")){
					for(i=0;i<6;i++) nodeinputs[i]=iomodel->Data("md.initialization.vz")[penta_vertex_ids[i]-1]*iomodel->Data("md.flowequation.borderFS")[penta_vertex_ids[i]-1];
					this->inputs->AddInput(new PentaInput(VzFSEnum,nodeinputs,P1Enum));
					for(i=0;i<6;i++) nodeinputs[i]=iomodel->Data("md.initialization.vz")[penta_vertex_ids[i]-1]*(1-iomodel->Data("md.flowequation.borderFS")[penta_vertex_ids[i]-1]);
					this->inputs->AddInput(new PentaInput(VzHOEnum,nodeinputs,P1Enum));
				}
				else{
					for(i=0;i<6;i++)nodeinputs[i]=0;
					this->inputs->AddInput(new PentaInput(VzFSEnum,nodeinputs,P1Enum));
					this->inputs->AddInput(new PentaInput(VzHOEnum,nodeinputs,P1Enum));
				}
			}
			if((IoCodeToEnumElementEquation(reCast<int>(iomodel->Data("md.flowequation.element_equation")[index])))==SSAFSApproximationEnum){
				/*Create VzSSA and VzFS Enums*/
				if(iomodel->Data("md.initialization.vz") && iomodel->Data("md.flowequation.borderFS")){
					for(i=0;i<6;i++) nodeinputs[i]=iomodel->Data("md.initialization.vz")[penta_vertex_ids[i]-1]*iomodel->Data("md.flowequation.borderFS")[penta_vertex_ids[i]-1];
					this->inputs->AddInput(new PentaInput(VzFSEnum,nodeinputs,P1Enum));
					for(i=0;i<6;i++) nodeinputs[i]=iomodel->Data("md.initialization.vz")[penta_vertex_ids[i]-1]*(1-iomodel->Data("md.flowequation.borderFS")[penta_vertex_ids[i]-1]);
					this->inputs->AddInput(new PentaInput(VzSSAEnum,nodeinputs,P1Enum));
				}
				else{
					for(i=0;i<6;i++)nodeinputs[i]=0;
					this->inputs->AddInput(new PentaInput(VzFSEnum,nodeinputs,P1Enum));
					this->inputs->AddInput(new PentaInput(VzSSAEnum,nodeinputs,P1Enum));
				}
			}
			break;
		default:
			/*No update for other solution types*/
			break;
	}
}
/*}}}*/
void       Penta::UpdateConstraintsExtrudeFromBase(void){/*{{{*/

	if(!IsOnBase()) return;

	int        extrusioninput;
	IssmDouble value,isonbase;

	this->parameters->FindParam(&extrusioninput,InputToExtrudeEnum);
	Input* input = inputs->GetInput(extrusioninput);      _assert_(extrusioninput);
	Input* onbase = inputs->GetInput(MeshVertexonbaseEnum); _assert_(onbase);

	GaussPenta* gauss=new GaussPenta();
	for(int iv=0;iv<this->NumberofNodes(this->element_type);iv++){
		gauss->GaussNode(this->element_type,iv);
		onbase->GetInputValue(&isonbase,gauss);
		if(isonbase==1.){
			input->GetInputValue(&value,gauss);
			this->nodes[iv]->ApplyConstraint(0,value);
		}
	}
	delete gauss;

}
/*}}}*/
void       Penta::UpdateConstraintsExtrudeFromTop(void){/*{{{*/

	if(!IsOnSurface()) return;

	int extrusioninput;
	int indices[3]={3,4,5};
	IssmDouble value;

	this->parameters->FindParam(&extrusioninput,InputToExtrudeEnum);
	Input* input = inputs->GetInput(extrusioninput); _assert_(extrusioninput);

	GaussPenta* gauss=new GaussPenta();
	for(int i=0;i<3;i++){
		gauss->GaussNode(P1Enum,indices[i]);
		input->GetInputValue(&value,gauss);
		this->nodes[indices[i]]->ApplyConstraint(0,value);
	}
	delete gauss;

}
/*}}}*/
int        Penta::UpdatePotentialUngrounding(IssmDouble* vertices_potentially_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf){/*{{{*/

	int i;
	int nflipped=0;

	/*Go through nodes, and whoever is on the potential_ungrounding, ends up in nodes_on_iceshelf: */
	for(i=0;i<NUMVERTICES;i++){
		if (reCast<bool,IssmDouble>(vertices_potentially_ungrounding[vertices[i]->Pid()])){
			vec_nodes_on_iceshelf->SetValue(vertices[i]->Pid(),-1.,INS_VAL);

			/*If node was not on ice shelf, we flipped*/
			if(nodes_on_iceshelf[vertices[i]->Pid()]>=0.){
				nflipped++;
			}
		}
	}
	return nflipped;
}
/*}}}*/
void       Penta::ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	PentaRef::GetInputDerivativeValue(dvalue,values,xyz_list,gauss,P1Enum);
}
/*}}}*/
void       Penta::ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss){/*{{{*/
	PentaRef::GetInputValue(pvalue,values,gauss,P1Enum);
}
/*}}}*/
int        Penta::VelocityInterpolation(void){/*{{{*/
	return PentaRef::VelocityInterpolation(this->element_type);
}
/*}}}*/
int        Penta::VertexConnectivity(int vertexindex){/*{{{*/
	_assert_(this->vertices);
	return this->vertices[vertexindex]->Connectivity();
}
/*}}}*/
void       Penta::VerticalSegmentIndices(int** pindices,int* pnumseg){/*{{{*/

	/*output*/
	int *indices = xNew<int>(3*2);
	indices[0*2 + 0] = 0; indices[0*2 + 1] = 3;
	indices[1*2 + 0] = 1; indices[1*2 + 1] = 4;
	indices[2*2 + 0] = 2; indices[2*2 + 1] = 5;

	/*Assign output pointers*/
	*pindices = indices;
	*pnumseg  = 3;
}
/*}}}*/
void       Penta::VerticalSegmentIndicesBase(int** pindices,int* pnumseg){/*{{{*/

	PentaRef::VerticalSegmentIndicesBase(pindices,pnumseg,this->GetElementType());

	return;
}
/*}}}*/
void       Penta::ViscousHeating(IssmDouble* pphi,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble phi;
	IssmDouble viscosity;
	IssmDouble epsilon[6];

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->StrainRateFS(&epsilon[0],xyz_list,(GaussPenta*)gauss,vx_input,vy_input,vz_input);
	this->material->ViscosityFS(&viscosity,3,xyz_list,(GaussPenta*)gauss,vx_input,vy_input,vz_input);
	GetPhi(&phi,&epsilon[0],viscosity);

	/*Assign output pointer*/
	*pphi = phi;
}
/*}}}*/
void       Penta::ZeroLevelsetCoordinates(IssmDouble** pxyz_zero,IssmDouble* xyz_list,int levelsetenum){/*{{{*/
	/*Compute portion of the element that is grounded*/ 

	int         normal_orientation=0;
	IssmDouble  s1,s2;
	IssmDouble  levelset[NUMVERTICES];
	IssmDouble* xyz_zero = xNew<IssmDouble>(4*3);

	/*Recover parameters and values*/
	GetInputListOnVertices(&levelset[0],levelsetenum);

	if(levelset[0]*levelset[1]>0.){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2
		/*Portion of the segments*/
		s1=levelset[2]/(levelset[2]-levelset[1]);
		s2=levelset[2]/(levelset[2]-levelset[0]);

		if(levelset[2]<0.) normal_orientation=1; //orientation of quadrangle at base and top, depending on distribution of levelsetfunction
		/*New point 1*/
		xyz_zero[3*normal_orientation+0]=xyz_list[2*3+0]+s1*(xyz_list[1*3+0]-xyz_list[2*3+0]);
		xyz_zero[3*normal_orientation+1]=xyz_list[2*3+1]+s1*(xyz_list[1*3+1]-xyz_list[2*3+1]);
		xyz_zero[3*normal_orientation+2]=xyz_list[2*3+2]+s1*(xyz_list[1*3+2]-xyz_list[2*3+2]);

		/*New point 0*/
		xyz_zero[3*(1-normal_orientation)+0]=xyz_list[2*3+0]+s2*(xyz_list[0*3+0]-xyz_list[2*3+0]);
		xyz_zero[3*(1-normal_orientation)+1]=xyz_list[2*3+1]+s2*(xyz_list[0*3+1]-xyz_list[2*3+1]);
		xyz_zero[3*(1-normal_orientation)+2]=xyz_list[2*3+2]+s2*(xyz_list[0*3+2]-xyz_list[2*3+2]);

		/*New point 3*/
		xyz_zero[3*(2+1-normal_orientation)+0]=xyz_list[5*3+0]+s1*(xyz_list[4*3+0]-xyz_list[5*3+0]);
		xyz_zero[3*(2+1-normal_orientation)+1]=xyz_list[5*3+1]+s1*(xyz_list[4*3+1]-xyz_list[5*3+1]);
		xyz_zero[3*(2+1-normal_orientation)+2]=xyz_list[5*3+2]+s1*(xyz_list[4*3+2]-xyz_list[5*3+2]);

		/*New point 4*/
		xyz_zero[3*(2+normal_orientation)+0]=xyz_list[5*3+0]+s2*(xyz_list[3*3+0]-xyz_list[5*3+0]);
		xyz_zero[3*(2+normal_orientation)+1]=xyz_list[5*3+1]+s2*(xyz_list[3*3+1]-xyz_list[5*3+1]);
		xyz_zero[3*(2+normal_orientation)+2]=xyz_list[5*3+2]+s2*(xyz_list[3*3+2]-xyz_list[5*3+2]);
	}
	else if(levelset[1]*levelset[2]>0.){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2
		/*Portion of the segments*/
		s1=levelset[0]/(levelset[0]-levelset[2]);
		s2=levelset[0]/(levelset[0]-levelset[1]);

		if(levelset[0]<0.) normal_orientation=1;
		/*New point 1*/
		xyz_zero[3*normal_orientation+0]=xyz_list[0*3+0]+s1*(xyz_list[2*3+0]-xyz_list[0*3+0]);
		xyz_zero[3*normal_orientation+1]=xyz_list[0*3+1]+s1*(xyz_list[2*3+1]-xyz_list[0*3+1]);
		xyz_zero[3*normal_orientation+2]=xyz_list[0*3+2]+s1*(xyz_list[2*3+2]-xyz_list[0*3+2]);

		/*New point 2*/
		xyz_zero[3*(1-normal_orientation)+0]=xyz_list[0*3+0]+s2*(xyz_list[1*3+0]-xyz_list[0*3+0]);
		xyz_zero[3*(1-normal_orientation)+1]=xyz_list[0*3+1]+s2*(xyz_list[1*3+1]-xyz_list[0*3+1]);
		xyz_zero[3*(1-normal_orientation)+2]=xyz_list[0*3+2]+s2*(xyz_list[1*3+2]-xyz_list[0*3+2]);

		/*New point 3*/
		xyz_zero[3*(2+1-normal_orientation)+0]=xyz_list[3*3+0]+s1*(xyz_list[5*3+0]-xyz_list[3*3+0]);
		xyz_zero[3*(2+1-normal_orientation)+1]=xyz_list[3*3+1]+s1*(xyz_list[5*3+1]-xyz_list[3*3+1]);
		xyz_zero[3*(2+1-normal_orientation)+2]=xyz_list[3*3+2]+s1*(xyz_list[5*3+2]-xyz_list[3*3+2]);

		/*New point 4*/
		xyz_zero[3*(2+normal_orientation)+0]=xyz_list[3*3+0]+s2*(xyz_list[4*3+0]-xyz_list[3*3+0]);
		xyz_zero[3*(2+normal_orientation)+1]=xyz_list[3*3+1]+s2*(xyz_list[4*3+1]-xyz_list[3*3+1]);
		xyz_zero[3*(2+normal_orientation)+2]=xyz_list[3*3+2]+s2*(xyz_list[4*3+2]-xyz_list[3*3+2]);
	}
	else if(levelset[0]*levelset[2]>0.){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2
		/*Portion of the segments*/
		s1=levelset[1]/(levelset[1]-levelset[0]);
		s2=levelset[1]/(levelset[1]-levelset[2]);

		if(levelset[1]<0.) normal_orientation=1;
		/*New point 0*/
		xyz_zero[3*normal_orientation+0]=xyz_list[1*3+0]+s1*(xyz_list[0*3+0]-xyz_list[1*3+0]);
		xyz_zero[3*normal_orientation+1]=xyz_list[1*3+1]+s1*(xyz_list[0*3+1]-xyz_list[1*3+1]);
		xyz_zero[3*normal_orientation+2]=xyz_list[1*3+2]+s1*(xyz_list[0*3+2]-xyz_list[1*3+2]);

		/*New point 2*/
		xyz_zero[3*(1-normal_orientation)+0]=xyz_list[1*3+0]+s2*(xyz_list[2*3+0]-xyz_list[1*3+0]);
		xyz_zero[3*(1-normal_orientation)+1]=xyz_list[1*3+1]+s2*(xyz_list[2*3+1]-xyz_list[1*3+1]);
		xyz_zero[3*(1-normal_orientation)+2]=xyz_list[1*3+2]+s2*(xyz_list[2*3+2]-xyz_list[1*3+2]);

		/*New point 3*/
		xyz_zero[3*(2+1-normal_orientation)+0]=xyz_list[4*3+0]+s1*(xyz_list[3*3+0]-xyz_list[4*3+0]);
		xyz_zero[3*(2+1-normal_orientation)+1]=xyz_list[4*3+1]+s1*(xyz_list[3*3+1]-xyz_list[4*3+1]);
		xyz_zero[3*(2+1-normal_orientation)+2]=xyz_list[4*3+2]+s1*(xyz_list[3*3+2]-xyz_list[4*3+2]);

		/*New point 4*/
		xyz_zero[3*(2+normal_orientation)+0]=xyz_list[4*3+0]+s2*(xyz_list[5*3+0]-xyz_list[4*3+0]);
		xyz_zero[3*(2+normal_orientation)+1]=xyz_list[4*3+1]+s2*(xyz_list[5*3+1]-xyz_list[4*3+1]);
		xyz_zero[3*(2+normal_orientation)+2]=xyz_list[4*3+2]+s2*(xyz_list[5*3+2]-xyz_list[4*3+2]);
	}
	else if(levelset[0]==0. && levelset[1]==0.){ //front is on point 0 and 1
		xyz_zero[3*0+0]=xyz_list[0*3+0];
		xyz_zero[3*0+1]=xyz_list[0*3+1];
		xyz_zero[3*0+2]=xyz_list[0*3+2];

		/*New point 2*/
		xyz_zero[3*1+0]=xyz_list[1*3+0];
		xyz_zero[3*1+1]=xyz_list[1*3+1];
		xyz_zero[3*1+2]=xyz_list[1*3+2];

		/*New point 3*/
		xyz_zero[3*2+0]=xyz_list[4*3+0];
		xyz_zero[3*2+1]=xyz_list[4*3+1];
		xyz_zero[3*2+2]=xyz_list[4*3+2];

		/*New point 4*/
		xyz_zero[3*3+0]=xyz_list[3*3+0];
		xyz_zero[3*3+1]=xyz_list[3*3+1];
		xyz_zero[3*3+2]=xyz_list[3*3+2];
	}
	else if(levelset[0]==0. && levelset[2]==0.){ //front is on point 0 and 1
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

		/*New point 4*/
		xyz_zero[3*3+0]=xyz_list[5*3+0];
		xyz_zero[3*3+1]=xyz_list[5*3+1];
		xyz_zero[3*3+2]=xyz_list[5*3+2];
	}
	else if(levelset[1]==0. && levelset[2]==0.){ //front is on point 0 and 1
		xyz_zero[3*0+0]=xyz_list[1*3+0];
		xyz_zero[3*0+1]=xyz_list[1*3+1];
		xyz_zero[3*0+2]=xyz_list[1*3+2];

		/*New point 2*/
		xyz_zero[3*1+0]=xyz_list[2*3+0];
		xyz_zero[3*1+1]=xyz_list[2*3+1];
		xyz_zero[3*1+2]=xyz_list[2*3+2];

		/*New point 3*/
		xyz_zero[3*2+0]=xyz_list[5*3+0];
		xyz_zero[3*2+1]=xyz_list[5*3+1];
		xyz_zero[3*2+2]=xyz_list[5*3+2];

		/*New point 4*/
		xyz_zero[3*3+0]=xyz_list[4*3+0];
		xyz_zero[3*3+1]=xyz_list[4*3+1];
		xyz_zero[3*3+2]=xyz_list[4*3+2];
	}
	else _error_("Case not covered");

	/*Assign output pointer*/
	*pxyz_zero= xyz_zero;
}
/*}}}*/

#ifdef _HAVE_GIAIVINS_
void       Penta::GiaDeflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt,IssmDouble* x,IssmDouble* y){/*{{{*/
	_error_("GIA deflection not implemented yet!");
}
/*}}}*/
#endif

#ifdef _HAVE_DAKOTA_
void       Penta::InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrows, int ncols, int name, int type){/*{{{*/

	int             i,t,row;
	IssmDouble      time;
	TransientInput *transientinput = NULL;
	IssmDouble      values[6];

	/*Check that name is an element input*/
	if(!IsInput(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	switch(type){

		case VertexEnum:
			/*Create transient input: */
			for(t=0;t<ncols;t++){ //ncols is the number of times

				/*create input values: */
				for(i=0;i<6;i++){
					row=this->vertices[i]->Sid();
					values[i]=matrix[ncols*row+t];
				}

				/*time:*/
				time=matrix[(nrows-1)*ncols+t];

				if(t==0) transientinput=new TransientInput(name);
				transientinput->AddTimeInput(new PentaInput(name,values,P1Enum),time);
				transientinput->Configure(parameters);
			}
			this->inputs->AddInput(transientinput);
			break;

		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}

}
/*}}}*/
void       Penta::InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){/*{{{*/

	int i,j;

	/*Check that name is an element input*/
	if(!IsInput(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	switch(type){

		case VertexEnum:

			/*New PentaInput*/
			IssmDouble values[6];

			/*Get values on the 6 vertices*/
			for (i=0;i<6;i++){
				values[i]=vector[this->vertices[i]->Sid()]; //careful, vector of values here is not parallel distributed, but serial distributed (from a serial Dakota core!)
			}

			/*Branch on the specified type of update: */
			switch(name){
				case ThicknessEnum:
					/*Update thickness + surface: assume bed is constant. On ice shelves, takes hydrostatic equilibrium*/
					IssmDouble  thickness[6];
					IssmDouble  thickness_init[6];
					IssmDouble  hydrostatic_ratio[6];
					IssmDouble  surface[6];
					IssmDouble  bed[6];

					/*retrieve inputs: */
					GetInputListOnVertices(&thickness_init[0],ThicknessEnum);
					GetInputListOnVertices(&hydrostatic_ratio[0],GeometryHydrostaticRatioEnum);
					GetInputListOnVertices(&bed[0],BaseEnum);
					GetInputListOnVertices(&surface[0],SurfaceEnum);

					/*build new thickness: */
//					for(j=0;j<6;j++)thickness[j]=values[j];

					/*build new bed and surface: */
					if (this->IsFloating()){
						/*hydrostatic equilibrium: */
						IssmDouble rho_ice,rho_water,di;
						rho_ice=this->matpar->GetMaterialParameter(MaterialsRhoIceEnum);
						rho_water=this->matpar->GetMaterialParameter(MaterialsRhoSeawaterEnum);

						di=rho_ice/rho_water;

						/*build new thickness: */
						for (j=0; j<6; j++) {
						/*  for observed/interpolated/hydrostatic thickness, remove scaling from any hydrostatic thickness  */
							if     (hydrostatic_ratio[j] >= 0.)
								thickness[j]=values[j]-(values[j]/thickness_init[j]-1.)*hydrostatic_ratio[j]*surface[j]/(1.-di);
						/*  for minimum thickness, don't scale  */
							else
								thickness[j]=thickness_init[j];

						/*  check the computed thickness and update bed  */
							if (thickness[j] < 0.)
								thickness[j]=1./(1.-di);
							bed[j]=surface[j]-thickness[j];
						}

//						for(j=0;j<6;j++){
//							surface[j]=(1-di)*thickness[j];
//							bed[j]=-di*thickness[j];
//						}
					}
					else{
						/*build new thickness: */
						for (j=0; j<6; j++) {
						/*  for observed thickness, use scaled value  */
							if(hydrostatic_ratio[j] >= 0.)
								thickness[j]=values[j];
						/*  for minimum thickness, don't scale  */
							else
								thickness[j]=thickness_init[j];
						}

						/*update bed on grounded ice: */
//						for(j=0;j<6;j++)surface[j]=bed[j]+thickness[j];
						for(j=0;j<6;j++)bed[j]=surface[j]-thickness[j];
					}

					/*Add new inputs: */
					this->inputs->AddInput(new PentaInput(ThicknessEnum,thickness,P1Enum));
					this->inputs->AddInput(new PentaInput(BaseEnum,bed,P1Enum));
					this->inputs->AddInput(new PentaInput(SurfaceEnum,surface,P1Enum));
					break;

				default:
					this->inputs->AddInput(new PentaInput(name,values,P1Enum));
			}
			break;

		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}

}
/*}}}*/
#endif
