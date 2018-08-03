/*!\file Tria.cpp
 * \brief: implementation of the Tria object
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
#include <math.h>
#include "../classes.h"
#include "../../shared/shared.h"
#ifdef _HAVE_GIAIVINS_
#include "../../modules/GiaDeflectionCorex/GiaDeflectionCorex.h"
#endif
/*}}}*/

/*Element macros*/
#define NUMVERTICES   3
#define NUMVERTICES1D 2

/*Constructors/destructor/copy*/
Tria::Tria(int tria_id, int tria_sid, int index, IoModel* iomodel,int nummodels)/*{{{*/
	:ElementHook(nummodels,index+1,NUMVERTICES,iomodel){

		/*id: */
		this->id  = tria_id;
		this->sid = tria_sid;

		//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
		this->parameters = NULL;

		/*intialize inputs: */
		this->inputs  = new Inputs();

		/*initialize pointers:*/
		this->nodes    = NULL;
		this->vertices = NULL;
		this->material = NULL;
		this->matpar   = NULL;
		if(nummodels>0){
			this->element_type_list=xNew<int>(nummodels);
			for(int i=0;i<nummodels;i++) this->element_type_list[i] = 0;
		}
		else this->element_type_list = NULL;

}
/*}}}*/
Tria::~Tria(){/*{{{*/
	this->parameters=NULL;
}
/*}}}*/
Object* Tria::copy() {/*{{{*/

	int i;
	Tria* tria=NULL;

	tria=new Tria();

	//deal with TriaRef mother class
	int nanalyses = this->numanalyses;
	if(nanalyses > 0){
		tria->element_type_list=xNew<int>(nanalyses);
		for(i=0;i<nanalyses;i++){
			if (this->element_type_list[i]) tria->element_type_list[i]=this->element_type_list[i];
			else tria->element_type_list[i] = 0;
		}
	}
	else tria->element_type_list = NULL;
	tria->element_type=this->element_type;
	tria->numanalyses=nanalyses;

	//deal with ElementHook mother class
	if (this->hnodes){
		tria->hnodes=xNew<Hook*>(tria->numanalyses);
		for(i=0;i<tria->numanalyses;i++){
			if (this->hnodes[i]) tria->hnodes[i] = (Hook*)(this->hnodes[i]->copy());
			else tria->hnodes[i] = NULL;
		}
	}
	else tria->hnodes = NULL;

	tria->hvertices = (Hook*)this->hvertices->copy();
	tria->hmaterial = (Hook*)this->hmaterial->copy();
	tria->hmatpar   = (Hook*)this->hmatpar->copy();
	tria->hneighbors = NULL;

	/*deal with Tria fields: */
	tria->id  = this->id;
	tria->sid = this->sid;
	if(this->inputs) tria->inputs = (Inputs*)(this->inputs->Copy());
	else tria->inputs=new Inputs();

	/*point parameters: */
	tria->parameters=this->parameters;

	/*recover objects: */
	if (this->nodes){
		unsigned int num_nodes = 3;
		tria->nodes = xNew<Node*>(num_nodes); //we cannot rely on an analysis_counter to tell us which analysis_type we are running, so we just copy the nodes.
		for(i=0;i<num_nodes;i++) if(this->nodes[i]) tria->nodes[i]=this->nodes[i]; else tria->nodes[i] = NULL;
	}
	else tria->nodes = NULL;
	
	tria->vertices = (Vertex**)this->hvertices->deliverp();
	tria->material = (Material*)this->hmaterial->delivers();
	tria->matpar   = (Matpar*)this->hmatpar->delivers();

	return tria;
}
/*}}}*/
void Tria::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/
	
	MARSHALLING_ENUM(TriaEnum);
	
	/*Call parent classes: */
	ElementHook::Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	Element::MarshallElement(pmarshalled_data,pmarshalled_data_size,marshall_direction,this->numanalyses);
	TriaRef::Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);

	vertices = (Vertex**)this->hvertices->deliverp();
	material = (Material*)this->hmaterial->delivers();
	matpar   = (Matpar*)this->hmatpar->delivers();

}
/*}}}*/

/*Other*/
void       Tria::AddBasalInput(int input_enum,IssmDouble* values, int interpolation_enum){/*{{{*/

	/*Call inputs method*/
	_assert_(this->inputs);
	
	int domaintype;
	parameters->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			this->inputs->AddInput(new TriaInput(input_enum,values,interpolation_enum));
			break;
		case Domain2DverticalEnum:{
			if(interpolation_enum==P1Enum){
				IssmDouble values2[NUMVERTICES]={0.};
				int        numindices;
				int       *indices = NULL;
				int        index = this->EdgeOnBaseIndex();
				NodeOnEdgeIndices(&numindices,&indices,index,this->FiniteElement());
				for(int i=0;i<numindices;i++){
					values2[indices[i]] = values[i];
				}
				this->inputs->AddInput(new TriaInput(input_enum,values2,interpolation_enum));
				xDelete<int>(indices);
			}
			else _error_("not implemented yet");
			}
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

}
/*}}}*/
void       Tria::AddInput(int input_enum,IssmDouble* values, int interpolation_enum){/*{{{*/

	/*Call inputs method*/
	_assert_(this->inputs);
	this->inputs->AddInput(new TriaInput(input_enum,values,interpolation_enum));
}
/*}}}*/
void       Tria::AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part){/*{{{*/

	bool       already = false;
	int        i,j;
	int        partition[NUMVERTICES];
	int        offsetsid[NUMVERTICES];
	int        offsetdof[NUMVERTICES];
	IssmDouble area;
	IssmDouble mean;

	/*First, get the area: */
	area=this->GetArea();

	/*Figure out the average for this element: */
	this->GetVerticesSidList(&offsetsid[0]);
	this->GetVertexPidList(&offsetdof[0]);
	mean=0;
	for(i=0;i<NUMVERTICES;i++){
		partition[i]=reCast<int>(qmu_part[offsetsid[i]]);
		mean=mean+1.0/NUMVERTICES*vertex_response[offsetdof[i]];
	}

	/*Add contribution: */
	for(i=0;i<NUMVERTICES;i++){
		already=false;
		for(j=0;j<i;j++){
			if (partition[i]==partition[j]){
				already=true;
				break;
			}
		}
		if(!already){
			partition_contributions->SetValue(partition[i],mean*area,ADD_VAL);
			partition_areas->SetValue(partition[i],area,ADD_VAL);
		};
	}
}
/*}}}*/
void       Tria::CalvingRateDev(){/*{{{*/

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

	/*Retrieve all inputs and parameters we will need*/
	Input* vx_input = inputs->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = inputs->GetInput(VyEnum); _assert_(vy_input);
	Input* gr_input = inputs->GetInput(MaskGroundediceLevelsetEnum); _assert_(gr_input);
	IssmDouble  B   = this->GetMaterialParameter(MaterialsRheologyBbarEnum);
	IssmDouble  n   = this->GetMaterialParameter(MaterialsRheologyNEnum);
	this->parameters->FindParam(&sigma_max_floating,CalvingStressThresholdFloatingiceEnum);
	this->parameters->FindParam(&sigma_max_grounded,CalvingStressThresholdGroundediceEnum);

	/* Start looping on the number of vertices: */
	GaussTria* gauss=new GaussTria();
	for(int iv=0;iv<NUMVERTICES;iv++){
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

		/*OLD (keep for a little bit)*/
		//sigma_max = 800.e+3; //IUGG previous test
		//sigma_max = 1000.e+3; //GRL
		//if(groundedice<0) sigma_max=150.e+3;

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
	this->inputs->AddInput(new TriaInput(CalvingratexEnum,&calvingratex[0],P1Enum));
	this->inputs->AddInput(new TriaInput(CalvingrateyEnum,&calvingratey[0],P1Enum));
	this->inputs->AddInput(new TriaInput(CalvingCalvingrateEnum,&calvingrate[0],P1Enum));

	/*Clean up and return*/
	delete gauss;
}
/*}}}*/
void       Tria::CalvingRateLevermann(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	GaussTria* gauss=NULL;
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
	Input* strainperpendicular_input=inputs->GetInput(StrainRateperpendicularEnum);					_assert_(strainperpendicular_input);
	Input* levermanncoeff_input=inputs->GetInput(CalvinglevermannCoeffEnum);                     _assert_(levermanncoeff_input);

	/* Start looping on the number of vertices: */
	gauss=new GaussTria();
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
	this->inputs->AddInput(new TriaInput(CalvingratexEnum,&calvingratex[0],P1Enum));
	this->inputs->AddInput(new TriaInput(CalvingrateyEnum,&calvingratey[0],P1Enum));
	this->inputs->AddInput(new TriaInput(CalvingCalvingrateEnum,&calvingrate[0],P1Enum));

	/*Clean up and return*/
	delete gauss;

}
/*}}}*/
IssmDouble Tria::CharacteristicLength(void){/*{{{*/

	return sqrt(2*this->GetArea());
}
/*}}}*/
void       Tria::ComputeBasalStress(Vector<IssmDouble>* eps){/*{{{*/
	_error_("Not Implemented yet");
}
/*}}}*/
void       Tria::ComputeDeviatoricStressTensor(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  viscosity;
	IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  tau_xx[NUMVERTICES];
	IssmDouble	tau_yy[NUMVERTICES];
	IssmDouble	tau_zz[NUMVERTICES]={0,0,0};
	IssmDouble  tau_xy[NUMVERTICES];
	IssmDouble	tau_xz[NUMVERTICES]={0,0,0};
	IssmDouble	tau_yz[NUMVERTICES]={0,0,0};
	IssmDouble  tau_e[NUMVERTICES];
	GaussTria*  gauss=NULL;
	int domaintype,dim=2;

	/*Get approximation*/
	int approximation;
	inputs->GetInputValue(&approximation,ApproximationEnum);

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	this->FindParam(&domaintype,DomainTypeEnum);
	Input* vx_input=inputs->GetInput(VxEnum);             _assert_(vx_input);
	Input* vy_input=inputs->GetInput(VyEnum);             _assert_(vy_input);

	/* Start looping on the number of vertices: */
	gauss=new GaussTria();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/*Compute strain rate and viscosity: */
		this->StrainRateSSA(&epsilon[0],&xyz_list[0][0],gauss,vx_input,vy_input);
		switch(approximation){
			case SSAApproximationEnum:
				this->material->ViscositySSA(&viscosity,dim,&xyz_list[0][0],gauss,vx_input,vy_input);
				break;
			case HOApproximationEnum:
				this->material->ViscosityHO(&viscosity,dim,&xyz_list[0][0],gauss,vx_input,vy_input);
				break;
			case FSApproximationEnum:
				this->material->ViscosityFS(&viscosity,dim,&xyz_list[0][0],gauss,vx_input,vy_input,NULL);
				break;
			default:
				_error_("not supported yet");
		}

		/*Compute Stress*/
		tau_xx[iv]=2*viscosity*epsilon[0]; // tau = nu eps
		tau_yy[iv]=2*viscosity*epsilon[1];
		tau_xy[iv]=2*viscosity*epsilon[2];
		tau_e[iv]=1/sqrt(2)*sqrt(pow(tau_xx[iv],2)+pow(tau_yy[iv],2)+2*pow(tau_xy[iv],2));
	}

	/*Add Stress tensor components into inputs*/
	this->inputs->AddInput(new TriaInput(DeviatoricStressxxEnum,&tau_xx[0],P1Enum));
	this->inputs->AddInput(new TriaInput(DeviatoricStressxyEnum,&tau_xy[0],P1Enum));
	this->inputs->AddInput(new TriaInput(DeviatoricStressxzEnum,&tau_xz[0],P1Enum));
	this->inputs->AddInput(new TriaInput(DeviatoricStressyyEnum,&tau_yy[0],P1Enum));
	this->inputs->AddInput(new TriaInput(DeviatoricStressyzEnum,&tau_yz[0],P1Enum));
	this->inputs->AddInput(new TriaInput(DeviatoricStresszzEnum,&tau_zz[0],P1Enum));
	this->inputs->AddInput(new TriaInput(DeviatoricStresseffectiveEnum,&tau_e[0],P1Enum));

	/*Clean up and return*/
	delete gauss;
}
/*}}}*/
void			Tria::ComputeEsaStrainAndVorticity(){ /*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  epsilon[4]; /* epsilon=[exx,eyy,exy+ (shear),exy- (rotation)];*/
	IssmDouble  strain_xx[NUMVERTICES];
	IssmDouble  strain_yy[NUMVERTICES];
	IssmDouble  strain_xy[NUMVERTICES];
	IssmDouble  vorticity_xy[NUMVERTICES];
	GaussTria*  gauss=NULL;
	
	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	Input* vx_input=this->GetInput(EsaEmotionEnum); _assert_(vx_input);
	Input* vy_input=this->GetInput(EsaNmotionEnum); _assert_(vy_input);
	
	/* Start looping on the number of vertices: */
	gauss=new GaussTria();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/*Compute strain rate and vorticity rate: */
		this->StrainRateESA(&epsilon[0],&xyz_list[0][0],gauss,vx_input,vy_input);

		/*Compute Stress*/
		strain_xx[iv]=epsilon[0];
		strain_yy[iv]=epsilon[1];
		strain_xy[iv]=epsilon[2];
		vorticity_xy[iv]=epsilon[3]; 
	}

	/*Add Stress tensor components into inputs*/
	this->inputs->AddInput(new TriaInput(EsaStrainratexxEnum,&strain_xx[0],P1Enum));
	this->inputs->AddInput(new TriaInput(EsaStrainrateyyEnum,&strain_yy[0],P1Enum));
	this->inputs->AddInput(new TriaInput(EsaStrainratexyEnum,&strain_xy[0],P1Enum));
	this->inputs->AddInput(new TriaInput(EsaRotationrateEnum,&vorticity_xy[0],P1Enum));

	/*Clean up and return*/
	delete gauss;
}
/*}}}*/
void       Tria::ComputeSigmaNN(){/*{{{*/

	if(!IsOnBase()){
		IssmDouble sigma_nn[3]={0.};
		this->inputs->AddInput(new TriaInput(SigmaNNEnum,&sigma_nn[0],P1Enum));
		return;
	}
	else{
		IssmDouble* xyz_list=NULL;
		IssmDouble *xyz_list_base=NULL;
		IssmDouble  pressure,viscosity;
		IssmDouble  sigma_nn[3];
		IssmDouble  sigma_xx,sigma_xy,sigma_yy;
		IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
		IssmDouble  base_normal[2]; 
		int domaintype,dim=2;

		/* Get node coordinates and dof list: */
		GetVerticesCoordinates(&xyz_list);
	   GetVerticesCoordinatesBase(&xyz_list_base);

		/*Retrieve all inputs we will be needing: */
		this->FindParam(&domaintype,DomainTypeEnum);
		if(domaintype==Domain2DhorizontalEnum) _error_("stress tensor calculation not supported for mesh of type " <<EnumToStringx(domaintype)<<", extrude mesh or call ComputeDeviatoricStressTensor");
		Input* pressure_input=inputs->GetInput(PressureEnum); _assert_(pressure_input);
		Input* vx_input=inputs->GetInput(VxEnum);             _assert_(vx_input);
		Input* vy_input=inputs->GetInput(VyEnum);             _assert_(vy_input);

		/* Start looping on the number of vertices: */
		Gauss* gauss = this->NewGauss();
		for(int i=0;i<NUMVERTICES;i++){
			gauss->GaussNode(P1Enum,i);

			/*Compute strain rate viscosity and pressure: */
			this->StrainRateSSA(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
			this->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,NULL);
			pressure_input->GetInputValue(&pressure,gauss);

			/*Compute Stress*/
			sigma_xx=2*viscosity*epsilon[0]-pressure; // sigma = nu eps - pressure
			sigma_yy=2*viscosity*epsilon[1]-pressure;
			sigma_xy=2*viscosity*epsilon[2];

			/*Get normal vector to the bed */
			NormalBase(&base_normal[0],xyz_list_base);

			/*Compute sigma_nn*/
			sigma_nn[i]=sigma_xx*base_normal[0]*base_normal[0] + 2*sigma_xy*base_normal[0]*base_normal[1] + sigma_yy*base_normal[1]*base_normal[1];
		}

		/*Add Stress tensor components into inputs*/
		this->inputs->AddInput(new TriaInput(SigmaNNEnum,&sigma_nn[0],P1Enum));

		/*Clean up and return*/
		xDelete<IssmDouble>(xyz_list);
		xDelete<IssmDouble>(xyz_list_base);
		delete gauss;
	}
}
/*}}}*/
void       Tria::ComputeStressTensor(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  pressure,viscosity;
	IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  sigma_xx[NUMVERTICES];
	IssmDouble	sigma_yy[NUMVERTICES];
	IssmDouble	sigma_zz[NUMVERTICES]={0,0,0};
	IssmDouble  sigma_xy[NUMVERTICES];
	IssmDouble	sigma_xz[NUMVERTICES]={0,0,0};
	IssmDouble	sigma_yz[NUMVERTICES]={0,0,0};
	GaussTria*  gauss=NULL;
	int domaintype,dim=2;

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	this->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype==Domain2DhorizontalEnum) _error_("stress tensor calculation not supported for mesh of type " <<EnumToStringx(domaintype)<<", extrude mesh or call ComputeDeviatoricStressTensor");
	Input* pressure_input=inputs->GetInput(PressureEnum); _assert_(pressure_input);
	Input* vx_input=inputs->GetInput(VxEnum);             _assert_(vx_input);
	Input* vy_input=inputs->GetInput(VyEnum);             _assert_(vy_input);

	/* Start looping on the number of vertices: */
	gauss=new GaussTria();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateSSA(&epsilon[0],&xyz_list[0][0],gauss,vx_input,vy_input);
		this->material->ViscositySSA(&viscosity,dim,&xyz_list[0][0],gauss,vx_input,vy_input);
		pressure_input->GetInputValue(&pressure,gauss);

		/*Compute Stress*/
		sigma_xx[iv]=2*viscosity*epsilon[0]-pressure; // sigma = nu eps - pressure
		sigma_yy[iv]=2*viscosity*epsilon[1]-pressure;
		sigma_xy[iv]=2*viscosity*epsilon[2];
	}

	/*Add Stress tensor components into inputs*/
	this->inputs->AddInput(new TriaInput(StressTensorxxEnum,&sigma_xx[0],P1Enum));
	this->inputs->AddInput(new TriaInput(StressTensorxyEnum,&sigma_xy[0],P1Enum));
	this->inputs->AddInput(new TriaInput(StressTensorxzEnum,&sigma_xz[0],P1Enum));
	this->inputs->AddInput(new TriaInput(StressTensoryyEnum,&sigma_yy[0],P1Enum));
	this->inputs->AddInput(new TriaInput(StressTensoryzEnum,&sigma_yz[0],P1Enum));
	this->inputs->AddInput(new TriaInput(StressTensorzzEnum,&sigma_zz[0],P1Enum));

	/*Clean up and return*/
	delete gauss;
}
/*}}}*/
void       Tria::Configure(Elements* elementsin, Loads* loadsin,Nodes* nodesin,Vertices *verticesin,Materials* materialsin, Parameters* parametersin){/*{{{*/

	/*go into parameters and get the analysis_counter: */
	int analysis_counter;
	parametersin->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Get Element type*/
	if (this->element_type_list) this->element_type=this->element_type_list[analysis_counter];

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective 
	 * datasets, using internal ids and offsets hidden in hooks: */
	if(this->hnodes){
		if (this->hnodes[analysis_counter]) this->hnodes[analysis_counter]->configure(nodesin);
		else this->hnodes[analysis_counter] = NULL;
	}
	else this->hnodes = NULL; 
	this->hvertices->configure(verticesin);
	this->hmaterial->configure(materialsin);
	this->hmatpar->configure(materialsin);

	/*Now, go pick up the objects inside the hooks: */
	if(this->hnodes && this->hnodes[analysis_counter]) this->nodes=(Node**)this->hnodes[analysis_counter]->deliverp();
	else this->nodes=NULL;
	this->vertices = (Vertex**)this->hvertices->deliverp();
	this->material = (Material*)this->hmaterial->delivers();
	this->matpar   = (Matpar*)this->hmatpar->delivers();

	/*point parameters to real dataset: */
	this->parameters=parametersin;

	/*get inputs configured too: */
	this->inputs->Configure(this->parameters);

}
/*}}}*/
void       Tria::ControlInputSetGradient(IssmDouble* gradient,int enum_type,int control_index){/*{{{*/

	int    vertexpidlist[NUMVERTICES];
	IssmDouble grad_list[NUMVERTICES];
	Input* grad_input=NULL;

	Input* input=inputs->GetInput(enum_type);
	if (!input) _error_("Input " << EnumToStringx(enum_type) << " not found");
	if (input->ObjectEnum()!=ControlInputEnum) _error_("Input " << EnumToStringx(enum_type) << " is not a ControlInput");

	GradientIndexing(&vertexpidlist[0],control_index);
	for(int i=0;i<NUMVERTICES;i++) grad_list[i]=gradient[vertexpidlist[i]];
	grad_input=new TriaInput(GradientEnum,grad_list,P1Enum);

	((ControlInput*)input)->SetGradient(grad_input);

}/*}}}*/
void       Tria::ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum){/*{{{*/

	Input* input=inputs->GetInput(control_enum);
	if (!input) _error_("Input " << EnumToStringx(control_enum) << " not found");
	if (input->ObjectEnum()!=ControlInputEnum) _error_("Input " << EnumToStringx(control_enum) << " is not a ControlInput");

	int         sidlist[NUMVERTICES];
	int         connectivity[NUMVERTICES];
	IssmPDouble values[NUMVERTICES];
	IssmPDouble gradients[NUMVERTICES]; 
	IssmDouble  value,gradient;

	this->GetVerticesConnectivityList(&connectivity[0]);
	this->GetVerticesSidList(&sidlist[0]);

	GaussTria* gauss=new GaussTria();
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
int        Tria::EdgeOnBaseIndex(void){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[3][2] = {{1,2},{2,0},{0,1}};

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);

	for(int i=0;i<3;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1.){
			return i;
		}
	}

	_printf_("list of vertices on bed: "<<values[0]<<" "<<values[1]<<" "<<values[2]);
	_error_("Could not find 2 vertices on bed");
}
/*}}}*/
void       Tria::EdgeOnBaseIndices(int* pindex1,int* pindex2){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[3][2] = {{1,2},{2,0},{0,1}};

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);

	for(int i=0;i<3;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1.){
			*pindex1 = indices[i][0];
			*pindex2 = indices[i][1];
			return;
		}
	}

	_printf_("list of vertices on bed: "<<values[0]<<" "<<values[1]<<" "<<values[2]);
	_error_("Could not find 2 vertices on bed");
}
/*}}}*/
int        Tria::EdgeOnSurfaceIndex(void){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[3][2] = {{1,2},{2,0},{0,1}};

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);

	for(int i=0;i<3;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1.){
			return i;
		}
	}

	_printf_("list of vertices on surface: "<<values[0]<<" "<<values[1]<<" "<<values[2]);
	_error_("Could not find 2 vertices on surface");
}
/*}}}*/
void       Tria::EdgeOnSurfaceIndices(int* pindex1,int* pindex2){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[3][2] = {{1,2},{2,0},{0,1}};

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);

	for(int i=0;i<3;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1.){
			*pindex1 = indices[i][0];
			*pindex2 = indices[i][1];
			return;
		}
	}

	_printf_("list of vertices on surface: "<<values[0]<<" "<<values[1]<<" "<<values[2]);
	_error_("Could not find 2 vertices on surface");
}
/*}}}*/
void       Tria::ElementResponse(IssmDouble* presponse,int response_enum){/*{{{*/

	switch(response_enum){
		case MaterialsRheologyBbarEnum:
			*presponse=this->material->GetBbar();
			break;

		case VelEnum:{

			/*Get input:*/
			IssmDouble vel;
			Input* vel_input;

			vel_input=this->inputs->GetInput(VelEnum); _assert_(vel_input);
			vel_input->GetInputAverage(&vel);

			/*Assign output pointers:*/
			*presponse=vel;}
			break;
		default:  
			_error_("Response type " << EnumToStringx(response_enum) << " not supported yet!");
	}

}
/*}}}*/
void       Tria::ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz){/*{{{*/

	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble xmin,ymin;
	IssmDouble xmax,ymax;

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	xmin=xyz_list[0][0]; xmax=xyz_list[0][0];
	ymin=xyz_list[0][1]; ymax=xyz_list[0][1];

	for(int i=1;i<NUMVERTICES;i++){
		if(xyz_list[i][0]<xmin) xmin=xyz_list[i][0];
		if(xyz_list[i][0]>xmax) xmax=xyz_list[i][0];
		if(xyz_list[i][1]<ymin) ymin=xyz_list[i][1];
		if(xyz_list[i][1]>ymax) ymax=xyz_list[i][1];
	}

	*hx=xmax-xmin;
	*hy=ymax-ymin;
	*hz=0.;
}
/*}}}*/
int        Tria::FiniteElement(void){/*{{{*/
	return this->element_type;
}
/*}}}*/
IssmDouble Tria::FloatingArea(void){/*{{{*/

	/*Intermediaries*/
	int         domaintype;
	IssmDouble  phi;
	IssmDouble *xyz_list  = NULL;

	if(!IsIceInElement())return 0.;

	/*Get problem dimension*/
	this->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain2DhorizontalEnum && domaintype!=Domain3DEnum) _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");

	this->GetVerticesCoordinates(&xyz_list);
	phi=this->GetGroundedPortion(xyz_list);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	return (1-phi)*this->GetArea();
}
/*}}}*/
void       Tria::FSContactMigration(Vector<IssmDouble>* vertexgrounded,Vector<IssmDouble>* vertexfloating){/*{{{*/

	if(!IsOnBase()) return;

	int approximation;
	inputs->GetInputValue(&approximation,ApproximationEnum);

	if(approximation==HOApproximationEnum || approximation==SSAApproximationEnum || approximation==SSAHOApproximationEnum){
		for(int i=0;i<NUMVERTICES;i++){
			vertexgrounded->SetValue(vertices[i]->Pid(),+9999.,INS_VAL);
			vertexfloating->SetValue(vertices[i]->Pid(),+9999.,INS_VAL);
		}
	}
	else{
		/*Intermediaries*/
		IssmDouble* xyz_list = NULL;
		IssmDouble* xyz_list_base = NULL;
		IssmDouble  pressure,water_pressure,sigma_nn,viscosity,bed,base;
		IssmDouble  bed_normal[2];
		IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
		IssmDouble  surface=0,value=0;
		bool grounded;

		/* Get node coordinates and dof list: */
		GetVerticesCoordinates(&xyz_list);
		GetVerticesCoordinatesBase(&xyz_list_base);

		/*Retrieve all inputs we will be needing: */
		Input* pressure_input = inputs->GetInput(PressureEnum); _assert_(pressure_input);
		Input* base_input     = inputs->GetInput(BaseEnum);     _assert_(base_input);
		Input* bed_input      = inputs->GetInput(BedEnum);      _assert_(bed_input);
		Input* vx_input       = inputs->GetInput(VxEnum);       _assert_(vx_input);
		Input* vy_input       = inputs->GetInput(VyEnum);       _assert_(vy_input);

		/*Create gauss point in the middle of the basal edge*/
		Gauss* gauss=NewGaussBase(1);
		gauss->GaussPoint(0);

		if(!IsFloating()){ 
			/*Check for basal force only if grounded and touching GL*/
			//		if(this->inputs->Min(MaskGroundediceLevelsetEnum)==0.){
			this->StrainRateSSA(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
			this->material->ViscosityFS(&viscosity,2,xyz_list,gauss,vx_input,vy_input,NULL);
			pressure_input->GetInputValue(&pressure, gauss);
			base_input->GetInputValue(&base, gauss); 

			/*Compute Stress*/
			IssmDouble sigma_xx=2.*viscosity*epsilon[0]-pressure;
			IssmDouble sigma_yy=2.*viscosity*epsilon[1]-pressure;
			IssmDouble sigma_xy=2.*viscosity*epsilon[2];

			/*Get normal vector to the bed */
			NormalBase(&bed_normal[0],xyz_list_base);

			/*basalforce*/
			sigma_nn = sigma_xx*bed_normal[0]*bed_normal[0] + sigma_yy*bed_normal[1]*bed_normal[1] + 2.*sigma_xy*bed_normal[0]*bed_normal[1];

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
			else         grounded=false;
		}
		for(int i=0;i<NUMVERTICES;i++){
			if(grounded) vertexgrounded->SetValue(vertices[i]->Pid(),+1.,INS_VAL);
			else         vertexfloating->SetValue(vertices[i]->Pid(),+1.,INS_VAL);
		}

		/*clean up*/
		delete gauss;
		xDelete<IssmDouble>(xyz_list);
		xDelete<IssmDouble>(xyz_list_base);
	}
}
/*}}}*/
IssmDouble Tria::GetArea(void){/*{{{*/

	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble x1,y1,x2,y2,x3,y3;

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x1=xyz_list[0][0]; y1=xyz_list[0][1];
	x2=xyz_list[1][0]; y2=xyz_list[1][1];
	x3=xyz_list[2][0]; y3=xyz_list[2][1];

	_assert_(x2*y3 - y2*x3 + x1*y2 - y1*x2 + x3*y1 - y3*x1>0);
	return (x2*y3 - y2*x3 + x1*y2 - y1*x2 + x3*y1 - y3*x1)/2;
}
/*}}}*/
IssmDouble Tria::GetArea3D(void){/*{{{*/

	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble x1,y1,z1,x2,y2,z2,x3,y3,z3;
	IssmDouble detm1,detm2,detm3;

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x1=xyz_list[0][0]; y1=xyz_list[0][1]; z1=xyz_list[0][2];
	x2=xyz_list[1][0]; y2=xyz_list[1][1]; z2=xyz_list[1][2];
	x3=xyz_list[2][0]; y3=xyz_list[2][1]; z3=xyz_list[2][2];

	detm1=x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2;
	detm2=y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2;
	detm3=x2*z1 - x1*z2 + x1*z3 - x3*z1 - x2*z3 + x3*z2;

	return sqrt(pow(detm1,2) + pow(detm2,2) + pow(detm3,2))/2;
}
/*}}}*/
IssmDouble Tria::GetAreaIce(void){/*{{{*/

	/*return area of element covered by ice*/
	/*Intermediaries*/
	int numiceverts;
	IssmDouble area_fraction;
	IssmDouble s[2]; // s:fraction of intersected triangle edges that lie inside ice
	int* indices=NULL;

	this->GetLevelsetIntersection(&indices, &numiceverts, s, MaskIceLevelsetEnum, 0.);

	switch (numiceverts){
		case 0: // no vertex has ice: element is ice free
			area_fraction=0.;
			break;
		case 1: // one vertex has ice: get area of triangle
			area_fraction=s[0]*s[1];
			break;
		case 2: // two vertices have ice: get area of quadrangle
			area_fraction=s[0]+s[1]-s[0]*s[1];
			break;
		case NUMVERTICES: // all vertices have ice: return triangle area
			area_fraction=1.;
			break;
		default:
			_error_("Wrong number of ice vertices in Tria::GetAreaIce!");
			break;
	}
	_assert_((area_fraction>=0.) && (area_fraction<=1.));

	xDelete<int>(indices);
	return area_fraction*this->GetArea();
}/*}}}*/
IssmDouble Tria::GetAreaSpherical(void){/*{{{*/

	bool spherical=true;
	IssmDouble llr_list[NUMVERTICES][3];
	IssmDouble x1,y1,z1,x2,y2,z2,x3,y3,z3;
	IssmDouble arc12,arc23,arc31,semi_peri,excess; 

	/*retrieve coordinates: lat,long,radius */
	::GetVerticesCoordinates(&llr_list[0][0],vertices,NUMVERTICES,spherical);
	x1=llr_list[0][0]/180*PI; y1=llr_list[0][1]/180*PI; z1=llr_list[0][2];
	x2=llr_list[1][0]/180*PI; y2=llr_list[1][1]/180*PI; z2=llr_list[1][2];
	x3=llr_list[2][0]/180*PI; y3=llr_list[2][1]/180*PI; z3=llr_list[2][2];

	/*compute great circle distance between vertices */
	arc12=2.*asin(sqrt(pow(sin((x2-x1)/2),2.0)+cos(x1)*cos(x2)*pow(sin((y2-y1)/2),2)));
	arc23=2.*asin(sqrt(pow(sin((x3-x2)/2),2.0)+cos(x2)*cos(x3)*pow(sin((y3-y2)/2),2)));
	arc31=2.*asin(sqrt(pow(sin((x1-x3)/2),2.0)+cos(x3)*cos(x1)*pow(sin((y1-y3)/2),2)));

	/*semi parameter */ 
	semi_peri=(arc12+arc23+arc31)/2; 

	/*spherical excess */
	excess=4.*atan(sqrt(tan(semi_peri/2)*tan((semi_peri-arc12)/2)*tan((semi_peri-arc23)/2)*tan((semi_peri-arc31)/2))); 

	/*area = excess*radius^2 */
	return excess*pow((z1+z2+z3)/3,2); 
}
/*}}}*/
void       Tria::GetAreaCoordinates(IssmDouble* area_coordinates,IssmDouble* xyz_zero,IssmDouble* xyz_list,int numpoints){/*{{{*/
	/*Computeportion of the element that is grounded*/ 

	int         i,j,k;
	IssmDouble  area_init,area_portion;
	IssmDouble  xyz_bis[NUMVERTICES][3];

	area_init=GetArea();

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
int        Tria::GetElementType(){/*{{{*/

	/*return TriaRef field*/
	return this->element_type;

}
/*}}}*/
void       Tria::GetGroundedPart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlyfloating){/*{{{*/
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
IssmDouble Tria::GetGroundedPortion(IssmDouble* xyz_list){/*{{{*/
	/*Computeportion of the element that is grounded*/ 

	bool              mainlyfloating = true;
	int               domaintype,index1,index2;
	const IssmPDouble epsilon        = 1.e-15;
	IssmDouble        phi,s1,s2,area_init,area_grounded;
	IssmDouble        gl[NUMVERTICES];
	IssmDouble        xyz_bis[3][3];

	/*Recover parameters and values*/
	parameters->FindParam(&domaintype,DomainTypeEnum);
	GetInputListOnVertices(&gl[0],MaskGroundediceLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	if(domaintype==Domain2DverticalEnum){
		this->EdgeOnBaseIndices(&index1,&index2);
		if(gl[index1]>0 && gl[index2]>0) phi=1; // All grounded
		else if(gl[index1]<0 && gl[index2]<0) phi=0; // All floating
		else if(gl[index1]<0 && gl[index2]>0){ //index2 grounded
			phi=1./(1.-gl[index1]/gl[index2]);
		}
		else if(gl[index2]<0 && gl[index1]>0){ //index1 grounded
			phi=1./(1.-gl[index2]/gl[index1]);
		}

	}
	else if(domaintype==Domain2DhorizontalEnum || domaintype==Domain3DEnum){
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
				xyz_bis[2][0]=*(xyz_list+3*2+0);
				xyz_bis[2][1]=*(xyz_list+3*2+1);
				xyz_bis[2][2]=*(xyz_list+3*2+2);

				/*Portion of the segments*/
				s1=gl[2]/(gl[2]-gl[1]);
				s2=gl[2]/(gl[2]-gl[0]);

				/*New point 1*/
				xyz_bis[1][0]=*(xyz_list+3*2+0)+s1*(*(xyz_list+3*1+0)-*(xyz_list+3*2+0));
				xyz_bis[1][1]=*(xyz_list+3*2+1)+s1*(*(xyz_list+3*1+1)-*(xyz_list+3*2+1));
				xyz_bis[1][2]=*(xyz_list+3*2+2)+s1*(*(xyz_list+3*1+2)-*(xyz_list+3*2+2));

				/*New point 0*/
				xyz_bis[0][0]=*(xyz_list+3*2+0)+s2*(*(xyz_list+3*0+0)-*(xyz_list+3*2+0));
				xyz_bis[0][1]=*(xyz_list+3*2+1)+s2*(*(xyz_list+3*0+1)-*(xyz_list+3*2+1));
				xyz_bis[0][2]=*(xyz_list+3*2+2)+s2*(*(xyz_list+3*0+2)-*(xyz_list+3*2+2));
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
			GetJacobianDeterminant(&area_init, xyz_list,NULL);
			GetJacobianDeterminant(&area_grounded, &xyz_bis[0][0],NULL);
			if(mainlyfloating==true) area_grounded=area_init-area_grounded;
			phi=area_grounded/area_init;
		}
	}
	else _error_("mesh type "<<EnumToStringx(domaintype)<<"not supported yet ");

	if(phi>1 || phi<0) _error_("Error. Problem with portion of grounded element: value should be between 0 and 1");

	return phi;
}
/*}}}*/
void       Tria::GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum){/*{{{*/
	
	/* Intermediaries */
	int i, dir,nrfrontnodes;
	IssmDouble  levelset[NUMVERTICES];
	int         indicesfront[NUMVERTICES];

	/*Recover parameters and values*/
	GetInputListOnVertices(&levelset[0],levelsetenum);

	/* Get nodes where there is no ice */
	nrfrontnodes=0;
	for(i=0;i<NUMVERTICES;i++){
		if(levelset[i]>=0.){
			indicesfront[nrfrontnodes]=i;
			nrfrontnodes++;
		}
	}

	_assert_(nrfrontnodes==2);

	/* arrange order of frontnodes such that they are oriented counterclockwise */
	if((NUMVERTICES+indicesfront[0]-indicesfront[1])%NUMVERTICES!=NUMVERTICES-1){
		int index=indicesfront[0];
		indicesfront[0]=indicesfront[1];
		indicesfront[1]=index;
	}	

	IssmDouble* xyz_front = xNew<IssmDouble>(3*nrfrontnodes);
	/* Return nodes */
	for(i=0;i<nrfrontnodes;i++){
		for(dir=0;dir<3;dir++){
			xyz_front[3*i+dir]=xyz_list[3*indicesfront[i]+dir];
		}
	}

	*pxyz_front=xyz_front;
}/*}}}*/
void       Tria::GetInputValue(IssmDouble* pvalue,Node* node,int enumtype){/*{{{*/

	Input* input=inputs->GetInput(enumtype);
	if(!input) _error_("No input of type " << EnumToStringx(enumtype) << " found in tria");

	GaussTria* gauss=new GaussTria();
	gauss->GaussVertex(this->GetNodeIndex(node));

	input->GetInputValue(pvalue,gauss);
	delete gauss;
}
/*}}}*/
void       Tria::GetLevelCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum,IssmDouble level){/*{{{*/

	/* Intermediaries */
	int i, dir,nrfrontnodes;
	IssmDouble  levelset[NUMVERTICES];
	int indicesfront[NUMVERTICES];

	/*Recover parameters and values*/
	GetInputListOnVertices(&levelset[0],levelsetenum);

	/* Get nodes where there is no ice */
	nrfrontnodes=0;
	for(i=0;i<NUMVERTICES;i++){
		if(levelset[i]==level){
			indicesfront[nrfrontnodes]=i;
			nrfrontnodes++;
		}
	}

	_assert_(nrfrontnodes==2);

	/* arrange order of frontnodes such that they are oriented counterclockwise */
	if((NUMVERTICES+indicesfront[0]-indicesfront[1])%NUMVERTICES!=NUMVERTICES-1){
		int index=indicesfront[0];
		indicesfront[0]=indicesfront[1];
		indicesfront[1]=index;
	}	

	IssmDouble* xyz_front = xNew<IssmDouble>(3*nrfrontnodes);
	/* Return nodes */
	for(i=0;i<nrfrontnodes;i++){
		for(dir=0;dir<3;dir++){
			xyz_front[3*i+dir]=xyz_list[3*indicesfront[i]+dir];
		}
	}

	*pxyz_front=xyz_front;

}/*}}}*/
void			Tria::GetLevelsetIntersection(int** pindices, int* pnumiceverts, IssmDouble* fraction, int levelset_enum, IssmDouble level){/*{{{*/
	
	/* GetLevelsetIntersection computes: 
	 * 1. indices of element, sorted in [iceverts, noiceverts] in counterclockwise fashion,
	 * 2. fraction of intersected triangle edges intersected by levelset, lying below level*/

	/*Intermediaries*/
	int i, numiceverts, numnoiceverts;
	int ind0, ind1, lastindex;
	int indices_ice[NUMVERTICES],indices_noice[NUMVERTICES];
	IssmDouble lsf[NUMVERTICES];
	int* indices = xNew<int>(NUMVERTICES);

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&lsf[0],levelset_enum);

	/* Determine distribution of ice over element.
	 * Exploit: ice/no-ice parts are connected, so find starting vertex of segment*/
	lastindex=0;
	for(i=0;i<NUMVERTICES;i++){ // go backwards along vertices, and check for sign change
		ind0=(NUMVERTICES-i)%NUMVERTICES;
		ind1=(ind0-1+NUMVERTICES)%NUMVERTICES;
		if((lsf[ind0]-level)*(lsf[ind1]-level)<=0.){ // levelset has been crossed, find last index belonging to segment
			if(lsf[ind1]==level) //if levelset intersects 2nd vertex, choose this vertex as last
				lastindex=ind1;
			else
				lastindex=ind0;
			break;
		}
	}

	numiceverts=0;
	numnoiceverts=0;
	for(i=0;i<NUMVERTICES;i++){
		ind0=(lastindex+i)%NUMVERTICES;
		if(lsf[i]<=level){
			indices_ice[numiceverts]=i;
			numiceverts++;
		}
		else{
			indices_noice[numnoiceverts]=i;
			numnoiceverts++;
		}
	}
	//merge indices 
	for(i=0;i<numiceverts;i++){indices[i]=indices_ice[i];}
	for(i=0;i<numnoiceverts;i++){indices[numiceverts+i]=indices_noice[i];}

	switch (numiceverts){
		case 0: // no vertex has ice: element is ice free, no intersection
			for(i=0;i<2;i++)
				fraction[i]=0.;
			break;
		case 1: // one vertex has ice:
			for(i=0;i<2;i++){
				fraction[i]=(level-lsf[indices[0]])/(lsf[indices[numiceverts+i]]-lsf[indices[0]]);
			}
			break;
		case 2: // two vertices have ice: fraction is computed from first ice vertex to last in CCW fashion
			for(i=0;i<2;i++){
				fraction[i]=(level-lsf[indices[i]])/(lsf[indices[numiceverts]]-lsf[indices[i]]);
			}
			break;
		case NUMVERTICES: // all vertices have ice: return triangle area
			for(i=0;i<2;i++)
				fraction[i]=1.;
			break;
		default:
			_error_("Wrong number of ice vertices in Tria::GetLevelsetIntersection!");
			break;
	}

	*pindices=indices;
	*pnumiceverts=numiceverts;
}
/*}}}*/
void       Tria::GetLevelsetPositivePart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlynegative,IssmDouble* gl){/*{{{*/
	
	/*Computeportion of the element that has a positive levelset*/ 

	bool               negative=true;
	int                point;
	const IssmPDouble  epsilon= 1.e-15;
	IssmDouble         f1,f2;

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	/*Check that not all nodes are positive or negative*/
	if(gl[0]>0 && gl[1]>0 && gl[2]>0){ // All positive
		point=0;
		f1=1.;
		f2=1.;
	}
	else if(gl[0]<0 && gl[1]<0 && gl[2]<0){ //All negative
		point=0;
		f1=0.;
		f2=0.;
	}
	else{
		if(gl[0]*gl[1]*gl[2]<0) negative=false;

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
	*mainlynegative=negative;
}
/*}}}*/
Node*      Tria::GetNode(int node_number){/*{{{*/
	_assert_(node_number>=0); 
	_assert_(node_number<this->NumberofNodes(this->element_type)); 
	return this->nodes[node_number];

}/*}}}*/
int        Tria::GetNodeIndex(Node* node){/*{{{*/

	_assert_(nodes);
	for(int i=0;i<NUMVERTICES;i++){
		if(node==nodes[i])
		 return i;
	}
	_error_("Node provided not found among element nodes");
}
/*}}}*/
int        Tria::GetNumberOfNodes(void){/*{{{*/
	if (this->nodes) return this->NumberofNodes(this->element_type);
	else return 0;
}
/*}}}*/
int        Tria::GetNumberOfNodes(int enum_type){/*{{{*/
	return this->NumberofNodes(enum_type);
}
/*}}}*/
int        Tria::GetNumberOfVertices(void){/*{{{*/
	return NUMVERTICES;
}
/*}}}*/
void       Tria::GetSolutionFromInputsOneDof(Vector<IssmDouble>* solution, int enum_type){/*{{{*/

	int        *doflist = NULL;
	IssmDouble  value;

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(this->element_type);

	/*Fetch dof list and allocate solution vector*/
	GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numnodes);

	/*Get inputs*/
	Input* enum_input=inputs->GetInput(enum_type); _assert_(enum_input);

	/*Ok, we have the values, fill in the array: */
	GaussTria* gauss=new GaussTria();
	for(int i=0;i<numnodes;i++){
		gauss->GaussNode(this->element_type,i);

		enum_input->GetInputValue(&value,gauss);
		values[i]=value;
	}

	solution->SetValues(numnodes,doflist,values,INS_VAL);

	/*Free ressources:*/
	xDelete<int>(doflist);
	xDelete<IssmDouble>(values);
	delete gauss;
}
/*}}}*/
void       Tria::GetVectorFromControlInputs(Vector<IssmDouble>* vector,int control_enum,int control_index,const char* data,bool onsid){/*{{{*/

	int vertexidlist[NUMVERTICES];
	Input *input=NULL;

	/*Get out if this is not an element input*/
	if(!IsInput(control_enum)) _error_("Enum "<<EnumToStringx(control_enum)<<" is not in IsInput");

	/*Prepare index list*/
	GradientIndexing(&vertexidlist[0],control_index,onsid);

	/*Get input (either in element or material)*/
	input=(Input*)this->inputs->GetInput(control_enum);   _assert_(input);

	/*Check that it is a ControlInput*/
	if (input->ObjectEnum()!=ControlInputEnum){
		_error_("input " << EnumToStringx(control_enum) << " is not a ControlInput");
	}

	((ControlInput*)input)->GetVectorFromInputs(vector,&vertexidlist[0],data);
}
/*}}}*/
void       Tria::GetVerticesCoordinatesBase(IssmDouble** pxyz_list){/*{{{*/

	int        indices[2];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Element XYZ list*/
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);

	/*Allocate Output*/
	IssmDouble* xyz_list_edge = xNew<IssmDouble>(2*3);
	this->EdgeOnBaseIndices(&indices[0],&indices[1]);
	for(int i=0;i<2;i++) for(int j=0;j<2;j++) xyz_list_edge[i*3+j]=xyz_list[indices[i]][j];

	/*Assign output pointer*/
	*pxyz_list = xyz_list_edge;

}/*}}}*/
void       Tria::GetVerticesCoordinatesTop(IssmDouble** pxyz_list){/*{{{*/

	int        indices[2];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Element XYZ list*/
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);

	/*Allocate Output*/
	IssmDouble* xyz_list_edge = xNew<IssmDouble>(2*3);
	this->EdgeOnSurfaceIndices(&indices[0],&indices[1]);
	for(int i=0;i<2;i++) for(int j=0;j<2;j++) xyz_list_edge[i*3+j]=xyz_list[indices[i]][j];

	/*Assign output pointer*/
	*pxyz_list = xyz_list_edge;

}/*}}}*/
IssmDouble Tria::GroundedArea(void){/*{{{*/

	/*Intermediaries*/
	int         domaintype;
	IssmDouble  phi;
	IssmDouble *xyz_list  = NULL;

	if(!IsIceInElement())return 0.;

	/*Get problem dimension*/
	this->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain2DhorizontalEnum && domaintype!=Domain3DEnum) _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");

	this->GetVerticesCoordinates(&xyz_list);
	phi=this->GetGroundedPortion(xyz_list);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	return phi*this->GetArea();
}
/*}}}*/
bool       Tria::HasEdgeOnBase(){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);
	sum = values[0]+values[1]+values[2];

	_assert_(sum==0. || sum==1. || sum==2.);

	if(sum==3.)  _error_("Two edges on bed not supported yet...");

	if(sum>1.){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
bool       Tria::HasEdgeOnSurface(){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);
	sum = values[0]+values[1]+values[2];

	_assert_(sum==0. || sum==1. || sum==2.);

	if(sum==3.)  _error_("Two edges on surface not supported yet...");

	if(sum>1.){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
IssmDouble Tria::IceMass(void){/*{{{*/

	IssmDouble rho_ice; 
	
	if(!IsIceInElement())return 0.; //do not contribute to the volume of the ice!

	/*recover ice density: */
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);

	return rho_ice*this->IceVolume();
}
/*}}}*/
IssmDouble Tria::IceVolume(void){/*{{{*/

	/*The volume of a truncated prism is area_base * 1/numedges sum(length of edges)*/

	/*Intermediaries*/
	int i, numiceverts;
	IssmDouble area_base,surface,base,Haverage;
	IssmDouble Haux[NUMVERTICES], surfaces[NUMVERTICES], bases[NUMVERTICES];
	IssmDouble s[2]; // s:fraction of intersected triangle edges, that lies inside ice
	int* indices=NULL;
	IssmDouble* H=NULL;

	if(!IsIceInElement())return 0.;

	int domaintype;
	parameters->FindParam(&domaintype,DomainTypeEnum);

	if(false && IsIcefront()){
		area_base=this->GetAreaIce();
		//Assumption: linear ice thickness profile on element. 
		//Hence ice thickness at intersection of levelset function with triangle edge is linear interpolation of ice thickness at vertices.
		this->GetLevelsetIntersection(&indices, &numiceverts, s, MaskIceLevelsetEnum, 0.);
		GetInputListOnVertices(&surfaces[0],SurfaceEnum);
		GetInputListOnVertices(&bases[0],BaseEnum);
		for(i=0;i<NUMVERTICES;i++) Haux[i]= surfaces[indices[i]]-bases[indices[i]]; //sort thicknesses in ice/noice
		int numthk=numiceverts+2;
		H=xNew<IssmDouble>(numthk);
		switch(numiceverts){
			case 1: // average over triangle 
				H[0]=Haux[0];
				H[1]=Haux[0]+s[0]*(Haux[1]-Haux[0]);
				H[2]=Haux[0]+s[1]*(Haux[2]-Haux[0]);
				break;
			case 2: // average over quadrangle
				H[0]=Haux[0];
				H[1]=Haux[1];
				H[2]=Haux[0]+s[0]*(Haux[2]-Haux[0]);
				H[3]=Haux[1]+s[1]*(Haux[2]-Haux[1]);
				break;
			default:
				_error_("Number of ice covered vertices wrong in Tria::IceVolume()");
				break;
		}
		Haverage=0.;
		for(i=0;i<numthk;i++)	Haverage+=H[i];
		Haverage/=IssmDouble(numthk);
	}
	else{
		/*First get back the area of the base*/
		area_base=this->GetArea();

		/*Now get the average height*/
		Input* surface_input = inputs->GetInput(SurfaceEnum); _assert_(surface_input);
		Input* base_input     = inputs->GetInput(BaseEnum);     _assert_(base_input);
		surface_input->GetInputAverage(&surface);
		base_input->GetInputAverage(&base);
		Haverage=surface-base;
	}

	/*Cleanup & return: */
	xDelete<int>(indices);
	xDelete<IssmDouble>(H);

	if(domaintype==Domain2DverticalEnum){
	  return area_base;
	}
	else{
	  return area_base*Haverage;
	}
}
/*}}}*/
IssmDouble Tria::IceVolumeAboveFloatation(void){/*{{{*/

	/*The volume above floatation: H + rho_water/rho_ice * bathymetry */
	IssmDouble rho_ice,rho_water;
	IssmDouble base,surface,bed,bathymetry;
	IssmDouble xyz_list[NUMVERTICES][3];

	if(!IsIceInElement() || IsFloating())return 0;

	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	rho_water=matpar->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Triangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2. * fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	/*Now get the average height and bathymetry*/
	Input* surface_input    = inputs->GetInput(SurfaceEnum);    _assert_(surface_input);
	Input* base_input        = inputs->GetInput(BaseEnum);        _assert_(base_input);
	Input* bed_input = inputs->GetInput(BedEnum); _assert_(bed_input);
	surface_input->GetInputAverage(&surface);
	base_input->GetInputAverage(&bed);
	bed_input->GetInputAverage(&bathymetry);
	
	/*Return: */
	return base*(surface-bed+min(rho_water/rho_ice*bathymetry,0.));
}
/*}}}*/
void       Tria::InputControlUpdate(IssmDouble scalar,bool save_parameter){/*{{{*/

	/*Intermediary*/
	int    num_controls;
	int*   control_type=NULL;
	Input* input=NULL;

	/*retrieve some parameters: */
	this->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	this->parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);

	for(int i=0;i<num_controls;i++){
		input=(Input*)this->inputs->GetInput(control_type[i]);   _assert_(input);
		if (input->ObjectEnum()!=ControlInputEnum){
			_error_("input " << EnumToStringx(control_type[i]) << " is not a ControlInput");
		}

		((ControlInput*)input)->UpdateValue(scalar);
		((ControlInput*)input)->Constrain();
		if (save_parameter) ((ControlInput*)input)->SaveValue();

	}

	/*Clean up and return*/
	xDelete<int>(control_type);
}
/*}}}*/
void       Tria::InputDepthAverageAtBase(int enum_type,int average_enum_type){/*{{{*/

	/*New input*/
	Input* oldinput=NULL;
	Input* newinput=NULL;

	/*copy input of enum_type*/
	oldinput=(Input*)this->inputs->GetInput(enum_type);
	if(!oldinput)_error_("could not find old input with enum: " << EnumToStringx(enum_type));
	newinput=(Input*)oldinput->copy();

	/*Assign new name (average)*/
	newinput->ChangeEnum(average_enum_type);

	/*Add new input to current element*/
	this->inputs->AddInput((Input*)newinput);
}
/*}}}*/
void       Tria::InputScale(int enum_type,IssmDouble scale_factor){/*{{{*/

	Input* input=NULL;

	/*Make a copy of the original input: */
	input=(Input*)this->inputs->GetInput(enum_type);
	if(!input)_error_("could not find old input with enum: " << EnumToStringx(enum_type));

	/*Scale: */
	input->Scale(scale_factor);
}
/*}}}*/
void       Tria::InputUpdateFromIoModel(int index, IoModel* iomodel){ //i is the element index/*{{{*/

	/*Intermediaries*/
	int        i,j;
	int        tria_vertex_ids[3];
	IssmDouble nodeinputs[3];
	IssmDouble cmmininputs[3];
	IssmDouble cmmaxinputs[3];
	bool       control_analysis   = false;
	int        num_control_type,num_responses;
	char**     controls = NULL;
	IssmDouble yts;

	/*Get parameters: */
	iomodel->FindConstant(&yts,"md.constants.yts"); 
	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");
	if(control_analysis) iomodel->FindConstant(&num_control_type,"md.inversion.num_control_parameters");
	if(control_analysis) iomodel->FindConstant(&num_responses,"md.inversion.num_cost_functions");

	/*Recover vertices ids needed to initialize inputs*/
	for(i=0;i<3;i++){ 
		tria_vertex_ids[i]=reCast<int>(iomodel->elements[3*index+i]); //ids for vertices are in the elements array from Matlab
	}

	/*Need to know the type of approximation for this element*/
	if(iomodel->Data("md.flowequation.element_equation")){
		this->inputs->AddInput(new IntInput(ApproximationEnum,IoCodeToEnumElementEquation(reCast<int>(iomodel->Data("md.flowequation.element_equation")[index]))));
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
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.balancethickness.thickening_rate")[tria_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i]/yts;
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i]/yts;
						this->inputs->AddInput(new ControlInput(BalancethicknessThickeningRateEnum,TriaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case VxEnum:
					if (iomodel->Data("md.initialization.vx")){
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.initialization.vx")[tria_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i]/yts;
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i]/yts;
						this->inputs->AddInput(new ControlInput(VxEnum,TriaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case VyEnum:
					if (iomodel->Data("md.initialization.vy")){
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.initialization.vy")[tria_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i]/yts;
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i]/yts;
						this->inputs->AddInput(new ControlInput(VyEnum,TriaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case ThicknessEnum:
					if(iomodel->Data("md.geometry.thickness")){
						for(j=0;j<NUMVERTICES;j++) nodeinputs[j]=iomodel->Data("md.geometry.thickness")[tria_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i];
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i];
						this->inputs->AddInput(new ControlInput(ThicknessEnum,TriaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case FrictionCoefficientEnum:
					if (iomodel->Data("md.friction.coefficient")){
						for(j=0;j<NUMVERTICES;j++)nodeinputs[j]=iomodel->Data("md.friction.coefficient")[tria_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i];
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i];
						this->inputs->AddInput(new ControlInput(FrictionCoefficientEnum,TriaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case MaterialsRheologyBbarEnum:
					if(iomodel->Data("md.materials.rheology_B")){
						for(j=0;j<NUMVERTICES;j++) nodeinputs[j]=iomodel->Data("md.materials.rheology_B")[tria_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i];
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i];
						this->inputs->AddInput(new ControlInput(MaterialsRheologyBbarEnum,TriaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				case DamageDbarEnum:
					if(iomodel->Data("md.damage.D")){
						for(j=0;j<NUMVERTICES;j++) nodeinputs[j]=iomodel->Data("md.damage.D")[tria_vertex_ids[j]-1];
						for(j=0;j<NUMVERTICES;j++)cmmininputs[j]=iomodel->Data("md.inversion.min_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i];
						for(j=0;j<NUMVERTICES;j++)cmmaxinputs[j]=iomodel->Data("md.inversion.max_parameters")[(tria_vertex_ids[j]-1)*num_control_type+i];
						this->inputs->AddInput(new ControlInput(DamageDbarEnum,TriaInputEnum,nodeinputs,cmmininputs,cmmaxinputs,i+1));
					}
					break;
				default:
					_error_("Control " << EnumToStringx(control) << " not implemented yet");
			}
		}
		for(i=0;i<num_control_type;i++) xDelete<char>(controls[i]);
		xDelete<char*>(controls);
	}

	/*DatasetInputs*/
	if (control_analysis && iomodel->Data("md.inversion.cost_functions_coefficients")){
	
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
			for(j=0;j<3;j++)nodeinputs[j]=iomodel->Data("md.inversion.cost_functions_coefficients")[(tria_vertex_ids[j]-1)*num_responses+i];
			datasetinput->AddInput(new TriaInput(InversionCostFunctionsCoefficientsEnum,nodeinputs,P1Enum),cost_functions_enums[i]);
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
void       Tria::InputUpdateFromSolutionOneDof(IssmDouble* solution,int enum_type){/*{{{*/

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
	this->inputs->AddInput(new TriaInput(enum_type,values,P1Enum));

	/*Free ressources:*/
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);
}
/*}}}*/
void       Tria::InputUpdateFromVector(IssmDouble* vector, int name, int type){/*{{{*/

	/*Check that name is an element input*/
	if(!IsInput(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	int         numnodes;
	int        *doflist = NULL;
	IssmDouble *values  = NULL;

	switch(type){
	case VertexPIdEnum: 
		values = xNew<IssmDouble>(NUMVERTICES);
		for(int i=0;i<NUMVERTICES;i++){
			values[i]=vector[this->vertices[i]->Pid()];
		}
		/*update input*/
		this->inputs->AddInput(new TriaInput(name,values,P1Enum));
		break;

	case VertexSIdEnum: 
		values = xNew<IssmDouble>(NUMVERTICES);
		for(int i=0;i<NUMVERTICES;i++){
			values[i]=vector[this->vertices[i]->Sid()];
		}
		/*update input*/
		this->inputs->AddInput(new TriaInput(name,values,P1Enum));
		break;

	case NodesEnum:
		/*Get number of nodes and dof list: */
		numnodes = this->NumberofNodes(this->element_type);
		values   = xNew<IssmDouble>(numnodes);
		GetDofList(&doflist,NoneApproximationEnum,GsetEnum);

		for(int i=0;i<numnodes;i++){
			values[i]=vector[doflist[i]];
			if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in vector");
			if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in vector");
		}
		this->inputs->AddInput(new TriaInput(name,values,this->element_type));
		break;

	case NodeSIdEnum:
		/*Get number of nodes and dof list: */
		numnodes = this->NumberofNodes(this->element_type);
		values   = xNew<IssmDouble>(numnodes);

		for(int i=0;i<numnodes;i++){
			values[i]=vector[nodes[i]->Sid()];
			if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in vector");
			if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in vector");
		}
		this->inputs->AddInput(new TriaInput(name,values,this->element_type));
		break;

	default:
		_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}

	/*Clean-up*/
	xDelete<int>(doflist);
	xDelete<IssmDouble>(values);

}
/*}}}*/
bool       Tria::IsFaceOnBoundary(void){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&values[0],MeshVertexonboundaryEnum);
	sum = values[0]+values[1]+values[2];

	_assert_(sum==0. || sum==1. || sum==2.);

	if(sum==3.)  _error_("Two edges on boundary not supported yet...");

	if(sum>1.){
		return true;
	}
	else{
		return false;
	}
}/*}}}*/
bool       Tria::IsIcefront(void){/*{{{*/

	bool isicefront;
	int i,nrice;
   IssmDouble ls[NUMVERTICES];

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&ls[0],MaskIceLevelsetEnum);

	/* If only one vertex has ice, there is an ice front here */
	isicefront=false;
	if(IsIceInElement()){
		nrice=0;       
		for(i=0;i<NUMVERTICES;i++)
			if(ls[i]<0.) nrice++;
		if(nrice==1) isicefront= true;
	}
	return isicefront;
}/*}}}*/
bool       Tria::IsNodeOnShelfFromFlags(IssmDouble* flags){/*{{{*/

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
bool       Tria::IsOnBase(){/*{{{*/

	int domaintype;
	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:
			return HasEdgeOnBase();
		case Domain2DhorizontalEnum:
			return true;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}
/*}}}*/
bool       Tria::IsOnSurface(){/*{{{*/

	int domaintype;
	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:
			return HasEdgeOnSurface();
		case Domain2DhorizontalEnum:
			return true;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}
/*}}}*/
bool       Tria::IsZeroLevelset(int levelset_enum){/*{{{*/

	bool iszerols;
	IssmDouble ls[NUMVERTICES];

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&ls[0],levelset_enum);

	/*If the level set is awlays <0, there is no ice front here*/
	iszerols= false;
	if(IsIceInElement()){
		if(ls[0]*ls[1]<0. || ls[0]*ls[2]<0. || (ls[0]*ls[1]*ls[2]==0. && ls[0]*ls[1]+ls[0]*ls[2]+ls[1]*ls[2]<=0.)){
			iszerols = true;
		}
	}

	return iszerols;
}
/*}}}*/
void       Tria::JacobianDeterminant(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetJacobianDeterminant(pJdet,xyz_list,(GaussTria*)gauss);

}
/*}}}*/
void       Tria::JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetSegmentJacobianDeterminant(pJdet,xyz_list_base,(GaussTria*)gauss);

}
/*}}}*/
void       Tria::JacobianDeterminantSurface(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetSegmentJacobianDeterminant(pJdet,xyz_list,(GaussTria*)gauss);

}
/*}}}*/
void       Tria::JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_top,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetSegmentJacobianDeterminant(pJdet,xyz_list_top,(GaussTria*)gauss);

}
/*}}}*/
IssmDouble Tria::Masscon(IssmDouble* levelset){ /*{{{*/


	/*intermediary: */
	IssmDouble* values=NULL;
	Input*      thickness_input=NULL;
	IssmDouble  thickness;
	IssmDouble  weight;
	IssmDouble  Jdet;
	IssmDouble  volume;
	IssmDouble  rho_ice;
	IssmDouble* xyz_list=NULL;
	int         point1;
	IssmDouble  fraction1,fraction2;
	bool        mainlynegative=true;
	
	/*Output:*/
	volume=0;

	/* Get node coordinates and dof list: */
	GetVerticesCoordinates(&xyz_list);

	/*Retrieve inputs required:*/
	thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);
	
	/*Retrieve material parameters: */
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);

	/*Retrieve values of the levelset defining the masscon: */
	values = xNew<IssmDouble>(NUMVERTICES);
	for(int i=0;i<NUMVERTICES;i++){
		values[i]=levelset[this->vertices[i]->Sid()];
	}
		
	/*Ok, use the level set values to figure out where we put our gaussian points:*/
	this->GetLevelsetPositivePart(&point1,&fraction1,&fraction2,&mainlynegative,values);
	Gauss* gauss = this->NewGauss(point1,fraction1,fraction2,mainlynegative,4);

	volume=0;

	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		this->JacobianDeterminant(&Jdet,xyz_list,gauss);
		thickness_input->GetInputValue(&thickness, gauss);

		volume+=thickness*gauss->weight*Jdet;
	}

	/* clean up and Return: */
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(values);
	delete gauss;
	return rho_ice*volume;
}
/*}}}*/
IssmDouble Tria::MassFlux(IssmDouble x1, IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id){/*{{{*/

	int        domaintype;
	IssmDouble mass_flux=0.;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble normal[2];
	IssmDouble length,rho_ice;
	IssmDouble h1,h2;
	IssmDouble vx1,vx2,vy1,vy2;
	GaussTria* gauss_1=NULL;
	GaussTria* gauss_2=NULL;

	/*Get material parameters :*/
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);

	/*First off, check that this segment belongs to this element: */
	if (segment_id!=this->id)_error_("error message: segment with id " << segment_id << " does not belong to element with id:" << this->id);

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*get area coordinates of 0 and 1 locations: */
	gauss_1=new GaussTria();
	gauss_1->GaussFromCoords(x1,y1,&xyz_list[0][0]);
	gauss_2=new GaussTria();
	gauss_2->GaussFromCoords(x2,y2,&xyz_list[0][0]);

	normal[0]=cos(atan2(x1-x2,y2-y1));
	normal[1]=sin(atan2(x1-x2,y2-y1));

	length=sqrt(pow(x2-x1,2)+pow(y2-y1,2));

	Input* thickness_input=inputs->GetInput(ThicknessEnum); _assert_(thickness_input);
	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	Input* vx_input=NULL;
	Input* vy_input=NULL;
	if(domaintype==Domain2DhorizontalEnum){
		vx_input=inputs->GetInput(VxEnum); _assert_(vx_input);
		vy_input=inputs->GetInput(VyEnum); _assert_(vy_input);
	}
	else{
		vx_input=inputs->GetInput(VxAverageEnum); _assert_(vx_input);
		vy_input=inputs->GetInput(VyAverageEnum); _assert_(vy_input);
	}

	thickness_input->GetInputValue(&h1, gauss_1);
	thickness_input->GetInputValue(&h2, gauss_2);
	vx_input->GetInputValue(&vx1,gauss_1);
	vx_input->GetInputValue(&vx2,gauss_2);
	vy_input->GetInputValue(&vy1,gauss_1);
	vy_input->GetInputValue(&vy2,gauss_2);

	mass_flux= rho_ice*length*(  
				(ONETHIRD*(h1-h2)*(vx1-vx2)+0.5*h2*(vx1-vx2)+0.5*(h1-h2)*vx2+h2*vx2)*normal[0]+
				(ONETHIRD*(h1-h2)*(vy1-vy2)+0.5*h2*(vy1-vy2)+0.5*(h1-h2)*vy2+h2*vy2)*normal[1]
				);

	/*clean up and return:*/
	delete gauss_1;
	delete gauss_2;
	return mass_flux;
}
/*}}}*/
IssmDouble Tria::MassFlux(IssmDouble* segment){/*{{{*/

	int        domaintype;
	IssmDouble mass_flux=0.;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble normal[2];
	IssmDouble length,rho_ice;
	IssmDouble x1,y1,x2,y2,h1,h2;
	IssmDouble vx1,vx2,vy1,vy2;
	GaussTria* gauss_1=NULL;
	GaussTria* gauss_2=NULL;

	/*Get material parameters :*/
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);

	/*First off, check that this segment belongs to this element: */
	if (reCast<int>(*(segment+4))!=this->id)_error_("error message: segment with id " << reCast<int>(*(segment+4)) << " does not belong to element with id:" << this->id);

	/*Recover segment node locations: */
	x1=*(segment+0); y1=*(segment+1); x2=*(segment+2); y2=*(segment+3);

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*get area coordinates of 0 and 1 locations: */
	gauss_1=new GaussTria();
	gauss_1->GaussFromCoords(x1,y1,&xyz_list[0][0]);
	gauss_2=new GaussTria();
	gauss_2->GaussFromCoords(x2,y2,&xyz_list[0][0]);

	normal[0]=cos(atan2(x1-x2,y2-y1));
	normal[1]=sin(atan2(x1-x2,y2-y1));

	length=sqrt(pow(x2-x1,2)+pow(y2-y1,2));

	Input* thickness_input=inputs->GetInput(ThicknessEnum); _assert_(thickness_input);
	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	Input* vx_input=NULL;
	Input* vy_input=NULL;
	if(domaintype==Domain2DhorizontalEnum){
		vx_input=inputs->GetInput(VxEnum); _assert_(vx_input);
		vy_input=inputs->GetInput(VyEnum); _assert_(vy_input);
	}
	else{
		vx_input=inputs->GetInput(VxAverageEnum); _assert_(vx_input);
		vy_input=inputs->GetInput(VyAverageEnum); _assert_(vy_input);
	}

	thickness_input->GetInputValue(&h1, gauss_1);
	thickness_input->GetInputValue(&h2, gauss_2);
	vx_input->GetInputValue(&vx1,gauss_1);
	vx_input->GetInputValue(&vx2,gauss_2);
	vy_input->GetInputValue(&vy1,gauss_1);
	vy_input->GetInputValue(&vy2,gauss_2);

	mass_flux= rho_ice*length*(  
				(ONETHIRD*(h1-h2)*(vx1-vx2)+0.5*h2*(vx1-vx2)+0.5*(h1-h2)*vx2+h2*vx2)*normal[0]+
				(ONETHIRD*(h1-h2)*(vy1-vy2)+0.5*h2*(vy1-vy2)+0.5*(h1-h2)*vy2+h2*vy2)*normal[1]
				);

	/*clean up and return:*/
	delete gauss_1;
	delete gauss_2;
	return mass_flux;
}
/*}}}*/
IssmDouble Tria::Misfit(int modelenum,int observationenum,int weightsenum){/*{{{*/

	/*Intermediaries*/
	IssmDouble model,observation,weight;
	IssmDouble Jdet;
	IssmDouble Jelem = 0;
	IssmDouble xyz_list[NUMVERTICES][3];
	GaussTria *gauss = NULL;

	/*If on water, return 0: */
	if(!IsIceInElement())return 0;

	/*Retrieve all inputs we will be needing: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* model_input=inputs->GetInput(modelenum);   _assert_(model_input);
	Input* observation_input=inputs->GetInput(observationenum);_assert_(observation_input);
	Input* weights_input     =inputs->GetInput(weightsenum);     _assert_(weights_input);

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		/* Get Jacobian determinant: */
		GetJacobianDeterminant(&Jdet, &xyz_list[0][0],gauss);

		/*Get parameters at gauss point*/
		model_input->GetInputValue(&model,gauss);
		observation_input->GetInputValue(&observation,gauss);
		weights_input->GetInputValue(&weight,gauss);

		/*compute misfit between model and observation */
		Jelem+=0.5*(model-observation)*(model-observation)*Jdet*weight*gauss->weight;
	}

	/* clean up and Return: */
	delete gauss;
	return Jelem;
}
/*}}}*/
IssmDouble Tria::MisfitArea(int weightsenum){/*{{{*/

	/*Intermediaries*/
	IssmDouble weight;
	IssmDouble Jdet;
	IssmDouble Jelem = 0;
	IssmDouble xyz_list[NUMVERTICES][3];
	GaussTria *gauss = NULL;

	/*If on water, return 0: */
	if(!IsIceInElement())return 0;

	/*Retrieve all inputs we will be needing: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* weights_input     =inputs->GetInput(weightsenum);     _assert_(weights_input);

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		/* Get Jacobian determinant: */
		GetJacobianDeterminant(&Jdet, &xyz_list[0][0],gauss);

		/*Get parameters at gauss point*/
		weights_input->GetInputValue(&weight,gauss);

		/*compute misfit between model and observation */
		Jelem+=Jdet*weight*gauss->weight;
	}

	/* clean up and Return: */
	delete gauss;
	return Jelem;
}
/*}}}*/
Gauss*     Tria::NewGauss(void){/*{{{*/
	return new GaussTria();
}
/*}}}*/
Gauss*     Tria::NewGauss(int order){/*{{{*/
	return new GaussTria(order);
}
/*}}}*/
Gauss*     Tria::NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order){/*{{{*/

	IssmDouble  area_coordinates[2][3];
	GetAreaCoordinates(&area_coordinates[0][0],xyz_list_front,xyz_list,2);
	return new GaussTria(area_coordinates,order);
}
/*}}}*/
Gauss*     Tria::NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order){/*{{{*/

	return new GaussTria(point1,fraction1,fraction2,mainlyfloating,order);
}
/*}}}*/
Gauss*     Tria::NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert){/*{{{*/

	IssmDouble  area_coordinates[2][3];
	GetAreaCoordinates(&area_coordinates[0][0],xyz_list_front,xyz_list,2);
	return new GaussTria(area_coordinates,order_vert);
}
/*}}}*/
Gauss*     Tria::NewGaussBase(int order){/*{{{*/

	int indices[2];
	this->EdgeOnBaseIndices(&indices[0],&indices[1]);
	return new GaussTria(indices[0],indices[1],order);
}
/*}}}*/
Gauss*     Tria::NewGaussTop(int order){/*{{{*/

	int indices[2];
	this->EdgeOnSurfaceIndices(&indices[0],&indices[1]);
	return new GaussTria(indices[0],indices[1],order);
}
/*}}}*/
void       Tria::NodalFunctions(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,this->element_type);

}
/*}}}*/
void       Tria::NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussTria*)gauss,this->element_type);

}
/*}}}*/
void       Tria::NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussTria*)gauss,this->VelocityInterpolation());

}
/*}}}*/
void       Tria::NodalFunctionsPressure(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,this->PressureInterpolation());

}
/*}}}*/
void       Tria::NodalFunctionsP1(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,P1Enum);

}
/*}}}*/
void       Tria::NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussTria*)gauss,P1Enum);

}
/*}}}*/
void       Tria::NodalFunctionsP2(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,P2Enum);

}
/*}}}*/
void       Tria::NodalFunctionsTensor(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,this->TensorInterpolation());

}
/*}}}*/
void       Tria::NodalFunctionsVelocity(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,this->VelocityInterpolation());

}
/*}}}*/
int        Tria::NodalValue(IssmDouble* pvalue, int index, int natureofdataenum){/*{{{*/

	int         found = 0;
	IssmDouble  value;
	Input      *data  = NULL;
	GaussTria  *gauss = NULL;

	/*First, serarch the input: */
	data=inputs->GetInput(natureofdataenum); 

	/*figure out if we have the vertex id: */
	found=0;
	for(int i=0;i<NUMVERTICES;i++){
		if(index==vertices[i]->Sid()){
			/*Do we have natureofdataenum in our inputs? :*/
			if(data){
				/*ok, we are good. retrieve value of input at vertex :*/
				gauss=new GaussTria(); gauss->GaussVertex(i);
				data->GetInputValue(&value,gauss);
				found=1;
				break;
			}
		}
	}

	/*clean-up*/
	delete gauss;

	if(found)*pvalue=value;
	return found;
}
/*}}}*/
void       Tria::NormalBase(IssmDouble* bed_normal,IssmDouble* xyz_list){/*{{{*/

	/*Build unit outward pointing vector*/
	IssmDouble vector[2];
	IssmDouble norm;

	vector[0]=xyz_list[1*3+0] - xyz_list[0*3+0];
	vector[1]=xyz_list[1*3+1] - xyz_list[0*3+1];

	norm=sqrt(vector[0]*vector[0] + vector[1]*vector[1]);

	bed_normal[0]= + vector[1]/norm;
	bed_normal[1]= - vector[0]/norm;
	_assert_(bed_normal[1]<0); 
}
/*}}}*/
void       Tria::NormalSection(IssmDouble* normal,IssmDouble* xyz_list){/*{{{*/

	/*Build unit outward pointing vector*/
	IssmDouble vector[2];
	IssmDouble norm;

	vector[0]=xyz_list[1*3+0] - xyz_list[0*3+0];
	vector[1]=xyz_list[1*3+1] - xyz_list[0*3+1];

	norm=sqrt(vector[0]*vector[0] + vector[1]*vector[1]);

	normal[0]= + vector[1]/norm;
	normal[1]= - vector[0]/norm;
}
/*}}}*/
void       Tria::NormalTop(IssmDouble* top_normal,IssmDouble* xyz_list){/*{{{*/

	/*Build unit outward pointing vector*/
	int index1,index2;
	IssmDouble vector[2];
	IssmDouble norm;

	this->EdgeOnSurfaceIndices(&index1,&index2);
	vector[0]=xyz_list[1*3+0] - xyz_list[0*3+0];
	vector[1]=xyz_list[1*3+1] - xyz_list[0*3+1];

	norm=sqrt(vector[0]*vector[0] + vector[1]*vector[1]);

	top_normal[0]= + vector[1]/norm;
	top_normal[1]= - vector[0]/norm;
	_assert_(top_normal[1]>0); 
}
/*}}}*/
int        Tria::ObjectEnum(void){/*{{{*/

	return TriaEnum;

}
/*}}}*/
int        Tria::NumberofNodesPressure(void){/*{{{*/
	return TriaRef::NumberofNodes(this->PressureInterpolation());
}
/*}}}*/
int        Tria::NumberofNodesVelocity(void){/*{{{*/
	return TriaRef::NumberofNodes(this->VelocityInterpolation());
}
/*}}}*/
void       Tria::PotentialUngrounding(Vector<IssmDouble>* potential_ungrounding){/*{{{*/

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

	/*go through vertices, and figure out which ones are grounded and want to unground: */
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
int        Tria::PressureInterpolation(void){/*{{{*/
	return TriaRef::PressureInterpolation(this->element_type);
}
/*}}}*/
void       Tria::ReduceMatrices(ElementMatrix* Ke,ElementVector* pe){/*{{{*/

	/*Static condensation if requested*/
	if(pe){
		if(this->element_type==MINIcondensedEnum){
			int indices[2]={6,7};
			pe->StaticCondensation(Ke,2,&indices[0]);
		}
		else if(this->element_type==P1bubblecondensedEnum){
			int size   = nodes[3]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int offset = 0;
			for(int i=0;i<3;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int* indices=xNew<int>(size);
			for(int i=0;i<size;i++) indices[i] = offset+i;
			pe->StaticCondensation(Ke,size,indices);
			xDelete<int>(indices);
		}
	}

	if(Ke){
		if(this->element_type==MINIcondensedEnum){
			int indices[2]={6,7};
			Ke->StaticCondensation(2,&indices[0]);
		}
		else if(this->element_type==P1bubblecondensedEnum){
			int size   = nodes[3]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int offset = 0;
			for(int i=0;i<3;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int* indices=xNew<int>(size);
			for(int i=0;i<size;i++) indices[i] = offset+i;
			Ke->StaticCondensation(size,indices);
			xDelete<int>(indices);
		}
	}


}
/*}}}*/
void       Tria::ResetFSBasalBoundaryCondition(void){/*{{{*/

	int numnodes = this->NumberofNodesVelocity();

	int          approximation;
	IssmDouble*  vertexonbase= NULL;
	IssmDouble   slope,groundedice;
	IssmDouble   xz_plane[6];

	/*For FS only: we want the CS to be tangential to the bedrock*/
	inputs->GetInputValue(&approximation,ApproximationEnum);
	if(!HasNodeOnBase() ||  approximation!=FSApproximationEnum) return;

	/*Get inputs*/
	Input* slope_input=inputs->GetInput(BedSlopeXEnum);                             _assert_(slope_input);
	Input* groundedicelevelset_input=inputs->GetInput(MaskGroundediceLevelsetEnum); _assert_(groundedicelevelset_input);
	vertexonbase = xNew<IssmDouble>(numnodes);
	this->GetInputListOnNodesVelocity(&vertexonbase[0],MeshVertexonbaseEnum);

	/*Loop over basal nodes and update their CS*/
	GaussTria* gauss = new GaussTria();
	for(int i=0;i<this->NumberofNodesVelocity();i++){

		if(vertexonbase[i]==1){
			gauss->GaussNode(this->VelocityInterpolation(),i);
			slope_input->GetInputValue(&slope,gauss);
			groundedicelevelset_input->GetInputValue(&groundedice,gauss);
			IssmDouble theta = atan(slope);

			/*New X axis                  New Z axis*/
			xz_plane[0]=cos(theta);       xz_plane[3]=0.;  
			xz_plane[1]=sin(theta);       xz_plane[4]=0.;  
			xz_plane[2]=0.;               xz_plane[5]=1.;          

			if(groundedice>=0){
				this->nodes[i]->DofInSSet(1); //vy
			}
			else{
				this->nodes[i]->DofInFSet(1); //vy
			}

			XZvectorsToCoordinateSystem(&this->nodes[i]->coord_system[0][0],&xz_plane[0]);
		}
	}

	/*cleanup*/
	xDelete<IssmDouble>(vertexonbase);
	delete gauss;
}
/*}}}*/
void       Tria::ResetHooks(){/*{{{*/

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
void       Tria::ResetLevelsetFromSegmentlist(IssmDouble* segments,int numsegments){/*{{{*/

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
	this->inputs->AddInput(new TriaInput(MaskIceLevelsetEnum,&ls[0],P1Enum));
}
/*}}}*/
void       Tria::SetClone(int* minranks){/*{{{*/

	_error_("not implemented yet");
}
/*}}}*/
void       Tria::SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index){/*{{{*/

	IssmDouble  values[NUMVERTICES];
	int         vertexpidlist[NUMVERTICES],control_init;


	/*Get Domain type*/
	int domaintype;
	parameters->FindParam(&domaintype,DomainTypeEnum);

	/*Specific case for depth averaged quantities*/
	control_init=control_enum;
	if(domaintype==Domain2DverticalEnum){
		if(control_enum==MaterialsRheologyBbarEnum){
			control_enum=MaterialsRheologyBEnum;
			if(!IsOnBase()) return;
		}
		if(control_enum==DamageDbarEnum){
			control_enum=DamageDEnum;
			if(!IsOnBase()) return;
		}
	}

	/*Get out if this is not an element input*/
	if(!IsInput(control_enum)) return;

	/*Prepare index list*/
	GradientIndexing(&vertexpidlist[0],control_index);

	/*Get values on vertices*/
	for(int i=0;i<NUMVERTICES;i++){
		values[i]=vector[vertexpidlist[i]];
	}
	Input* new_input = new TriaInput(control_enum,values,P1Enum);
	Input* input     = (Input*)this->inputs->GetInput(control_enum);   _assert_(input);
	if(input->ObjectEnum()!=ControlInputEnum){
		_error_("input " << EnumToStringx(control_enum) << " is not a ControlInput");
	}

	((ControlInput*)input)->SetInput(new_input);
}
/*}}}*/
void       Tria::SetCurrentConfiguration(Elements* elementsin, Loads* loadsin, Nodes* nodesin, Materials* materialsin, Parameters* parametersin){/*{{{*/

	/*go into parameters and get the analysis_counter: */
	int analysis_counter;
	parametersin->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Get Element type*/
	if(this->element_type_list) this->element_type=this->element_type_list[analysis_counter];

	/*Pick up nodes*/
	if(this->hnodes && this->hnodes[analysis_counter]){
		this->nodes=(Node**)this->hnodes[analysis_counter]->deliverp();
	}

}
/*}}}*/
Element*   Tria::SpawnBasalElement(void){/*{{{*/

	int index1,index2;
	int domaintype;

	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			return this;
		case Domain2DverticalEnum:
			_assert_(HasEdgeOnBase());
			this->EdgeOnBaseIndices(&index1,&index2);
			return SpawnSeg(index1,index2);
		default:
			_error_("not implemented yet");
	}
}
/*}}}*/
Seg*       Tria::SpawnSeg(int index1,int index2){/*{{{*/

	int analysis_counter;

	/*go into parameters and get the analysis_counter: */
	this->parameters->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Create Seg*/
	Seg* seg=new Seg();
	seg->id=this->id;
	seg->inputs=(Inputs*)this->inputs->SpawnSegInputs(index1,index2);
	seg->parameters=this->parameters;
	seg->element_type=P1Enum; //Only P1 CG for now (TO BE CHANGED)
	this->SpawnSegHook(xDynamicCast<ElementHook*>(seg),index1,index2);

	/*Spawn material*/
	seg->material=(Material*)this->material->copy2(seg);

	/*recover nodes, material and matpar: */
	seg->nodes    = (Node**)seg->hnodes[analysis_counter]->deliverp();
	seg->vertices = (Vertex**)seg->hvertices->deliverp();
	seg->matpar   = (Matpar*)seg->hmatpar->delivers();

	/*Return new Seg*/
	return seg;
}
/*}}}*/
Element*   Tria::SpawnTopElement(void){/*{{{*/

	int index1,index2;
	int domaintype;

	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			return this;
		case Domain2DverticalEnum:
			_assert_(HasEdgeOnSurface());
			this->EdgeOnSurfaceIndices(&index1,&index2);
			return SpawnSeg(index2,index1); //reverse order
		default:
			_error_("not implemented yet");
	}
}
/*}}}*/
void       Tria::StrainRateparallel(){/*{{{*/

	IssmDouble *xyz_list = NULL;
	IssmDouble  epsilon[3];
	GaussTria* gauss=NULL;
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

	/* Start looping on the number of vertices: */
	gauss=new GaussTria();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/* Get the value we need*/
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vel=vx*vx+vy*vy;

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateSSA(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
		strainxx=epsilon[0];
		strainyy=epsilon[1];
		strainxy=epsilon[2];

		/*strainparallel= Strain rate along the ice flow direction */
		strainparallel[iv]=(vx*vx*(strainxx)+vy*vy*(strainyy)+2*vy*vx*strainxy)/(vel+1.e-14);
	}

	/*Add input*/
	this->inputs->AddInput(new TriaInput(StrainRateparallelEnum,&strainparallel[0],P1Enum));

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
}
/*}}}*/
void       Tria::StrainRateperpendicular(){/*{{{*/

	IssmDouble *xyz_list = NULL;
	GaussTria* gauss=NULL;
	IssmDouble  epsilon[3];
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

	/* Start looping on the number of vertices: */
	gauss=new GaussTria();
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/* Get the value we need*/
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vel=vx*vx+vy*vy;

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateSSA(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
		strainxx=epsilon[0];
		strainyy=epsilon[1];
		strainxy=epsilon[2];

		/*strainperpendicular= Strain rate perpendicular to the ice flow direction */
		strainperpendicular[iv]=(vx*vx*(strainyy)+vy*vy*(strainxx)-2*vy*vx*strainxy)/(vel+1.e-14);
	}

	/*Add input*/
	this->inputs->AddInput(new TriaInput(StrainRateperpendicularEnum,&strainperpendicular[0],P1Enum));

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
}
/*}}}*/
IssmDouble Tria::SurfaceArea(void){/*{{{*/

	IssmDouble S;
	IssmDouble normal[3];
	IssmDouble v13[3],v23[3];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*If on water, return 0: */
	if(!IsIceInElement()) return 0.;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	for(int i=0;i<3;i++){
		v13[i]=xyz_list[0][i]-xyz_list[2][i];
		v23[i]=xyz_list[1][i]-xyz_list[2][i];
	}

	normal[0]=v13[1]*v23[2]-v13[2]*v23[1];
	normal[1]=v13[2]*v23[0]-v13[0]*v23[2];
	normal[2]=v13[0]*v23[1]-v13[1]*v23[0];

	S = 0.5 * sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

	/*Return: */
	return S;
}
/*}}}*/
int        Tria::TensorInterpolation(void){/*{{{*/
	return TriaRef::TensorInterpolation(this->element_type);
}
/*}}}*/
IssmDouble Tria::TimeAdapt(void){/*{{{*/

	/*intermediary: */
	int    i;
	IssmDouble C,dt;
	IssmDouble dx,dy;
	IssmDouble maxx,minx;
	IssmDouble maxy,miny;
	IssmDouble maxabsvx,maxabsvy;
	IssmDouble xyz_list[NUMVERTICES][3];

	/*get CFL coefficient:*/
	this->parameters->FindParam(&C,TimesteppingCflCoefficientEnum);

	/*Get for Vx and Vy, the max of abs value: */
	maxabsvx = this->inputs->MaxAbs(VxEnum);
	maxabsvy = this->inputs->MaxAbs(VyEnum);

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	minx=xyz_list[0][0];
	maxx=xyz_list[0][0];
	miny=xyz_list[0][1];
	maxy=xyz_list[0][1];

	for(i=1;i<NUMVERTICES;i++){
		if (xyz_list[i][0]<minx)minx=xyz_list[i][0];
		if (xyz_list[i][0]>maxx)maxx=xyz_list[i][0];
		if (xyz_list[i][1]<miny)miny=xyz_list[i][1];
		if (xyz_list[i][1]>maxy)maxy=xyz_list[i][1];
	}
	dx=maxx-minx;
	dy=maxy-miny;

	/*CFL criterion: */
	dt=C/(maxabsvx/dx+maxabsvy/dy);

	return dt;
}
/*}}}*/
IssmDouble Tria::TotalFloatingBmb(void){/*{{{*/

	/*The fbmb[kg yr-1] of one element is area[m2] * melting_rate [kg m^-2 yr^-1]*/
	int        point1;
	bool       mainlyfloating;
	IssmDouble fbmb=0;
	IssmDouble rho_ice,fraction1,fraction2,floatingmelt,Jdet;
	IssmDouble Total_Fbmb=0;
	IssmDouble xyz_list[NUMVERTICES][3];
	Gauss*     gauss     = NULL;

   if(!IsIceInElement())return 0;

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
		this->JacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		floatingmelt_input->GetInputValue(&floatingmelt,gauss);
		fbmb+=floatingmelt*Jdet*gauss->weight;
	}

   Total_Fbmb=rho_ice*fbmb;	        // from volume to mass

	/*Return: */
	delete gauss;
	return Total_Fbmb;
}
/*}}}*/
IssmDouble Tria::TotalGroundedBmb(void){/*{{{*/

	/*The gbmb[kg yr-1] of one element is area[m2] * gounded melting rate [kg m^-2 yr^-1]*/
	int        point1;
	bool       mainlyfloating;
	IssmDouble gbmb=0;
	IssmDouble rho_ice,fraction1,fraction2,groundedmelt,Jdet;
	IssmDouble Total_Gbmb=0;
	IssmDouble xyz_list[NUMVERTICES][3];
	Gauss*     gauss     = NULL;

   if(!IsIceInElement())return 0;

	/*Get material parameters :*/
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	Input* groundedmelt_input = this->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(groundedmelt_input); 
	Input* gllevelset_input = this->GetInput(MaskGroundediceLevelsetEnum); _assert_(gllevelset_input);
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	this->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating);
	/* Start  looping on the number of gaussian points: */
	gauss = this->NewGauss(point1,fraction1,fraction2,mainlyfloating,2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);
		this->JacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		groundedmelt_input->GetInputValue(&groundedmelt,gauss);
		gbmb+=groundedmelt*Jdet*gauss->weight;
	}

   Total_Gbmb=rho_ice*gbmb;	        // from volume to mass

	/*Return: */
	delete gauss;
	return Total_Gbmb;
}
/*}}}*/
IssmDouble Tria::TotalSmb(void){/*{{{*/

	/*The smb[kg yr-1] of one element is area[m2] * smb [kg m^-2 yr^-1]*/
	IssmDouble base,smb,rho_ice;
	IssmDouble Total_Smb=0;
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Get material parameters :*/
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);

   if(!IsIceInElement())return 0;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Triangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2. * fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));	// area of element in m2

	/*Now get the average SMB over the element*/
	Input* smb_input = inputs->GetInput(SmbMassBalanceEnum); _assert_(smb_input);
	smb_input->GetInputAverage(&smb);																								// average smb on element in m ice s-1
   Total_Smb=rho_ice*base*smb;																											// smb on element in kg s-1

	/*Return: */
	return Total_Smb;
}
/*}}}*/
void       Tria::Update(int index, IoModel* iomodel,int analysis_counter,int analysis_type,int finiteelement_type){/*{{{*/

	/*Intermediaries*/
	int  numnodes;
	int* tria_node_ids = NULL;

	/*Checks if debuging*/
	_assert_(iomodel->elements);

	/*Recover element type*/
	this->element_type_list[analysis_counter]=finiteelement_type;

	/*Recover nodes ids needed to initialize the node hook.*/
	switch(finiteelement_type){
		case P1Enum:
			numnodes        = 3;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->nodecounter+iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->nodecounter+iomodel->elements[3*index+2];
			break;
		case P1DGEnum:
			numnodes        = 3;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+3*index+1;
			tria_node_ids[1]=iomodel->nodecounter+3*index+2;
			tria_node_ids[2]=iomodel->nodecounter+3*index+3;
			break;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			numnodes        = 4;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->nodecounter+iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->nodecounter+iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->nodecounter+iomodel->numberofvertices+index+1;
			break;
		case P2Enum:
			numnodes        = 6;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->nodecounter+iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->nodecounter+iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;
			break;
		case P2bubbleEnum: case P2bubblecondensedEnum:
			numnodes        = 7;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->nodecounter+iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->nodecounter+iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;
			tria_node_ids[6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+index+1;
			break;
		case P1P1Enum: case P1P1GLSEnum:
			numnodes        = 6;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->nodecounter+iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->nodecounter+iomodel->elements[3*index+2];

			tria_node_ids[3]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elements[3*index+0];
			tria_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elements[3*index+1];
			tria_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elements[3*index+2];
			break;
		case MINIEnum: case MINIcondensedEnum:
			numnodes       = 7;
			tria_node_ids  = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->nodecounter+iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->nodecounter+iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->nodecounter+iomodel->numberofvertices+index+1;

			tria_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[3*index+0];
			tria_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[3*index+1];
			tria_node_ids[6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[3*index+2];
			break;
		case TaylorHoodEnum:
		case XTaylorHoodEnum:
			numnodes        = 9;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->nodecounter+iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->nodecounter+iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;

			tria_node_ids[6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[3*index+0];
			tria_node_ids[7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[3*index+1];
			tria_node_ids[8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[3*index+2];
			break;
		case LATaylorHoodEnum:
			numnodes        = 6;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->nodecounter+iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->nodecounter+iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;
			break;
		case CrouzeixRaviartEnum:
			numnodes        = 10;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->nodecounter+iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->nodecounter+iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;
			tria_node_ids[6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+index+1;

			tria_node_ids[7]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements+3*index+1;
			tria_node_ids[8]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements+3*index+2;
			tria_node_ids[9]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements+3*index+3;
			break;
		case LACrouzeixRaviartEnum:
			numnodes        = 7;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->nodecounter+iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->nodecounter+iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->nodecounter+iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;
			tria_node_ids[6]=iomodel->nodecounter+iomodel->numberofvertices+iomodel->numberofedges+index+1;
			break;
		default:
			_error_("Finite element "<<EnumToStringx(finiteelement_type)<<" not supported yet");
	}

	/*hooks: */
	this->SetHookNodes(tria_node_ids,numnodes,analysis_counter); this->nodes=NULL;
	xDelete<int>(tria_node_ids);

	/*Fill with IoModel*/
	this->InputUpdateFromIoModel(index,iomodel);
}
/*}}}*/
void       Tria::UpdateConstraintsExtrudeFromBase(void){/*{{{*/

	if(!HasEdgeOnBase()) return;

	int        extrusioninput;
	IssmDouble value,isonbase;

	this->parameters->FindParam(&extrusioninput,InputToExtrudeEnum);
	Input* input = inputs->GetInput(extrusioninput);      _assert_(input);
	Input* onbase = inputs->GetInput(MeshVertexonbaseEnum); _assert_(onbase);

	GaussTria* gauss=new GaussTria();
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
void       Tria::UpdateConstraintsExtrudeFromTop(void){/*{{{*/

	if(!HasEdgeOnSurface()) return;

	int extrusioninput;
	int indices[2];
	IssmDouble value;

	this->parameters->FindParam(&extrusioninput,InputToExtrudeEnum);
	Input* input = inputs->GetInput(extrusioninput); _assert_(input);
	this->EdgeOnSurfaceIndices(&indices[0],&indices[1]);

	GaussTria* gauss=new GaussTria();
	for(int i=0;i<2;i++){
		gauss->GaussNode(P1Enum,indices[i]);
		input->GetInputValue(&value,gauss);
		this->nodes[indices[i]]->ApplyConstraint(0,value);
	}
	delete gauss;

}
/*}}}*/
int        Tria::UpdatePotentialUngrounding(IssmDouble* vertices_potentially_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf){/*{{{*/

	int i;
	int nflipped=0;

	/*Go through nodes, and whoever is on the potential_ungrounding, ends up in nodes_on_iceshelf: */
	for(i=0;i<3;i++){
		if (reCast<bool>(vertices_potentially_ungrounding[vertices[i]->Pid()])){
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
void       Tria::ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	TriaRef::GetInputDerivativeValue(dvalue,values,xyz_list,gauss,P1Enum);
}
/*}}}*/
void       Tria::ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss){/*{{{*/
	TriaRef::GetInputValue(pvalue,values,gauss,P1Enum);
}
/*}}}*/
int        Tria::VelocityInterpolation(void){/*{{{*/
	return TriaRef::VelocityInterpolation(this->element_type);
}
/*}}}*/
int        Tria::VertexConnectivity(int vertexindex){/*{{{*/
	_assert_(this->vertices);
	return this->vertices[vertexindex]->Connectivity();
}
/*}}}*/
void       Tria::WriteLevelsetSegment(DataSet* segments){/*{{{*/

	if(!this->IsZeroLevelset(MaskIceLevelsetEnum)) return;

	IssmDouble* xyz_list_zero = NULL;
	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);
	this->ZeroLevelsetCoordinates(&xyz_list_zero,&xyz_list[0][0], MaskIceLevelsetEnum);
	if(xyz_list_zero){
		IssmDouble x[2];
		IssmDouble y[2];
		x[0] = xyz_list_zero[0*3 + 0]; x[1] = xyz_list_zero[1*3 + 0];
		y[0] = xyz_list_zero[0*3 + 1]; y[1] = xyz_list_zero[1*3 + 1];
		segments->AddObject(new Contour<IssmDouble>(segments->Size()+1,2,&x[0],&y[0],false));
	}
	xDelete<IssmDouble>(xyz_list_zero);

//	IssmDouble ls[NUMVERTICES];
//	IssmDouble  xyz_list[NUMVERTICES][3];
//
//	if(IsIceInElement()){
//
//		/*Retrieve all inputs and parameters*/
//		GetInputListOnVertices(&ls[0],levelset_enum);
//
//		/*If the level set is awlays <0, there is no ice front here*/
//		bool iszerols= false;
//		if(IsIceInElement()){
//			if(ls[0]*ls[1]<0. || ls[0]*ls[2]<0. || (ls[0]*ls[1]*ls[2]==0. && ls[0]*ls[1]+ls[0]*ls[2]+ls[1]*ls[2]<=0.)){
//				iszerols = true;
//			}
//		}
//
//		if(iszerols){
//			/*OK we have one segment!*/
//			::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
//			int count = 0;
//			IssmDouble x[2];
//			IssmDouble y[2];
//
//			for(int i=0;i<NUMVERTICES,i++){
//				int index1 = i;
//				int index1 = (i+1)%3;
//				if(ls[index1]<=0 && ls[index2]>=0){
//
//				}
//
//			}
//			Contour* segment = new Contour<IssmDouble>(segment->Size()+1,2,x,y,false);
//		}
//
//	}
//
//	_error_("STOP");
}
/*}}}*/
void       Tria::ZeroLevelsetCoordinates(IssmDouble** pxyz_zero,IssmDouble* xyz_list,int levelsetenum){/*{{{*/
	/* Return coordinates where levelset intersects element edges.
	 * Attention: In case that no intersection exists, NULL pointer is returned.*/

	/*Intermediaries*/
	const int dim=3;
	int numiceverts;
	int i, n, e, counter;
	IssmDouble s[2];
	int* indices=NULL;
	IssmDouble* xyz_zero=NULL;

	this->GetLevelsetIntersection(&indices, &numiceverts, s,levelsetenum,0.);
	
	//TODO: check if for 2 iceverts front segment is oriented in CCW way
	
	if(numiceverts>0) xyz_zero=xNew<IssmDouble>(2*dim);
	if((numiceverts>0)&&(numiceverts<NUMVERTICES)){
		counter=0;
		for(i=0;i<numiceverts;i++){	// iterate over ice vertices
			for(n=numiceverts;n<NUMVERTICES;n++){ // iterate over no-ice vertices
				for(e=0;e<dim;e++){ // spatial direction
					int ind_ice		=dim*indices[i]+e;
					int ind_noice	=dim*indices[n]+e;
					int ind			=dim*counter+e;
					xyz_zero[ind]=xyz_list[ind_ice]+s[counter]*(xyz_list[ind_noice]-xyz_list[ind_ice]);
				}
				counter++;
			}
		}
	}
	else if(numiceverts==NUMVERTICES){ //NUMVERTICES ice vertices: calving front lies on element edge
		IssmDouble lsf[NUMVERTICES];
		this->GetInputListOnVertices(&lsf[0],levelsetenum);
		counter=0;
		for(i=0;i<NUMVERTICES;i++){
			if(lsf[indices[i]]==0.){
				for(e=0;e<dim;e++)	xyz_zero[dim*counter+e]=xyz_list[dim*indices[i]+e];
				counter++;
			}
			if(counter==2) break;
		}
	}
	_assert_(counter==2);

	/*Cleanup & return*/
	xDelete<int>(indices);
	*pxyz_zero=xyz_zero;
}
/*}}}*/

#ifdef _HAVE_GIAIVINS_
void       Tria::GiaDeflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt,IssmDouble* x, IssmDouble* y){/*{{{*/

	int i;
	int gsize;
	IssmDouble xi,yi,ri,re,area;
	IssmDouble x0,y0;
	IssmDouble xyz_list[NUMVERTICES][3];

	/*thickness averages: */
	IssmDouble* hes=NULL;
	IssmDouble* times=NULL;
	IssmDouble  currenttime;
	int         numtimes;
	Input* thickness_input=NULL;

	/*gia solution parameters:*/
	int cross_section_shape=0;

	/*gia material parameters: */
	IssmDouble lithosphere_shear_modulus;
	IssmDouble lithosphere_density;
	IssmDouble mantle_shear_modulus;
	IssmDouble mantle_density;
	Input* mantle_viscosity_input=NULL;
	IssmDouble mantle_viscosity;
	Input* lithosphere_thickness_input=NULL;
	IssmDouble lithosphere_thickness;

	/*ice properties: */
	IssmDouble rho_ice;

	/*constants: */
	IssmDouble yts;

	/*output: */
	IssmDouble  wi;
	IssmDouble  dwidt;

	/*arguments to GiaDeflectionCorex: */
	GiaDeflectionCoreArgs arguments;

	/*how many dofs are we working with here? */
	this->parameters->FindParam(&gsize,MeshNumberofverticesEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);

	/*recover gia solution parameters: */
	this->parameters->FindParam(&cross_section_shape,GiaCrossSectionShapeEnum);

	/*what time is it? :*/
	this->parameters->FindParam(&currenttime,TimeEnum);

	/*recover material parameters: */
	lithosphere_shear_modulus=matpar->GetMaterialParameter(MaterialsLithosphereShearModulusEnum);
	lithosphere_density=matpar->GetMaterialParameter(MaterialsLithosphereDensityEnum);
	mantle_shear_modulus=matpar->GetMaterialParameter(MaterialsMantleShearModulusEnum);
	mantle_density=matpar->GetMaterialParameter(MaterialsMantleDensityEnum);
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);

	/*pull thickness averages: */
	thickness_input=inputs->GetInput(ThicknessEnum); 
	if (!thickness_input)_error_("thickness input needed to compute gia deflection!");
	thickness_input->GetInputUpToCurrentTimeAverages(&hes,&times,&numtimes,currenttime);

	/*recover mantle viscosity: */
	mantle_viscosity_input=inputs->GetInput(GiaMantleViscosityEnum);
	if (!mantle_viscosity_input)_error_("mantle viscosity input needed to compute gia deflection!");
	mantle_viscosity_input->GetInputAverage(&mantle_viscosity);

	/*recover lithosphere thickness: */
	lithosphere_thickness_input=inputs->GetInput(GiaLithosphereThicknessEnum);
	if (!lithosphere_thickness_input)_error_("lithosphere thickness input needed to compute gia deflection!");
	lithosphere_thickness_input->GetInputAverage(&lithosphere_thickness);

	/*pull area of this Tria: */
	area=this->GetArea();

	/*element radius: */
	re=sqrt(area/PI);

	/*figure out gravity center of our element: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x0=(xyz_list[0][0]+xyz_list[1][0]+xyz_list[2][0])/3.0;
	y0=(xyz_list[0][1]+xyz_list[1][1]+xyz_list[2][1])/3.0;

	/*start loading GiaDeflectionCore arguments: */
	arguments.re=re;
	arguments.hes=hes;
	arguments.times=times;
	arguments.numtimes=numtimes;
	arguments.currenttime=currenttime;
	arguments.lithosphere_shear_modulus=lithosphere_shear_modulus;
	arguments.lithosphere_density=lithosphere_density;
	arguments.mantle_shear_modulus=mantle_shear_modulus;
	arguments.mantle_viscosity=mantle_viscosity;
	arguments.mantle_density=mantle_density;
	arguments.lithosphere_thickness=lithosphere_thickness;
	arguments.rho_ice=rho_ice;
	arguments.idisk=this->id;
	arguments.iedge=cross_section_shape;
	arguments.yts=yts;

	for(i=0;i<gsize;i++){
		/*compute distance from the center of the tria to the vertex i: */
		xi=x[i]; yi=y[i];
		ri=sqrt(pow(xi-x0,2)+pow(yi-y0,2));

		/*load ri onto arguments for this vertex i: */
		arguments.ri=ri;

		/*for this Tria, compute contribution to rebound at vertex i: */
		GiaDeflectionCorex(&wi,&dwidt,&arguments);

		/*plug value into solution vector: */
		wg->SetValue(i,wi,ADD_VAL);
		dwgdt->SetValue(i,dwidt,ADD_VAL);

	}

	/*Free ressources: */
	xDelete<IssmDouble>(hes);
	xDelete<IssmDouble>(times);

	return;
}
/*}}}*/
#endif
#ifdef _HAVE_ESA_
void    Tria::EsaGeodetic2D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* xx,IssmDouble* yy){ /*{{{*/

	/*diverse:*/
	int gsize;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble area;
	IssmDouble earth_radius = 6371012.0;	// Earth's radius [m]
	IssmDouble I;		//ice/water loading 
	IssmDouble rho_ice, rho_earth;

	/*precomputed elastic green functions:*/
	IssmDouble* U_elastic_precomputed = NULL;
	IssmDouble* H_elastic_precomputed = NULL;
	int         M;
	
	/*computation of Green functions:*/
	IssmDouble* U_elastic= NULL;
	IssmDouble* N_elastic= NULL;
	IssmDouble* E_elastic= NULL;
	
	/*optimization:*/
	bool store_green_functions=false;

	/*Compute ice thickness change: */
	Input*	deltathickness_input=inputs->GetInput(EsaDeltathicknessEnum); 
	if (!deltathickness_input)_error_("delta thickness input needed to compute elastic adjustment!");
	deltathickness_input->GetInputAverage(&I);
		
	/*early return if we are not on the (ice) loading point: */ 
	if(I==0) return; 

	/*recover material parameters: */
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	rho_earth=matpar->GetMaterialParameter(MaterialsEarthDensityEnum);

	/*how many dofs are we working with here? */
	this->parameters->FindParam(&gsize,MeshNumberofverticesEnum);

	/*compute area of element:*/
	area=GetArea();

	/*figure out gravity center of our element (Cartesian): */
	IssmDouble x_element, y_element; 
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x_element=(xyz_list[0][0]+xyz_list[1][0]+xyz_list[2][0])/3.0;
	y_element=(xyz_list[0][1]+xyz_list[1][1]+xyz_list[2][1])/3.0;

	/*recover elastic Green's functions for displacement:*/
	DoubleVecParam* U_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(EsaUElasticEnum)); _assert_(U_parameter);
	DoubleVecParam* H_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(EsaHElasticEnum)); _assert_(H_parameter);
	U_parameter->GetParameterValueByPointer(&U_elastic_precomputed,&M);
	H_parameter->GetParameterValueByPointer(&H_elastic_precomputed,&M);

	/*initialize: */
	U_elastic=xNewZeroInit<IssmDouble>(gsize);
	N_elastic=xNewZeroInit<IssmDouble>(gsize);
	E_elastic=xNewZeroInit<IssmDouble>(gsize);

	int* indices=xNew<int>(gsize);
	IssmDouble* U_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* N_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* E_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble dx, dy; 
	IssmDouble dist, alpha, ang;
	IssmDouble N_azim, E_azim;

	for(int i=0;i<gsize;i++){

		indices[i]=i; 
		IssmDouble N_azim=0; 
		IssmDouble E_azim=0;

		/*Compute alpha angle between centroid and current vertex: */
		dx = x_element - xx[i];		dy = y_element - yy[i]; 
		dist = sqrt(pow(dx,2)+pow(dy,2));						// distance between vertex and elemental centroid [m] 
		alpha = dist*360.0/(2*PI*earth_radius) * PI/180.0;	// [in radians] 360 degree = 2*pi*earth_radius 

		/*Compute azimuths, both north and east components: */
		ang = PI/2 - atan2(dy,dx);		// this is bearing angle! 
		N_azim = cos(ang); 
		E_azim = sin(ang); 

		/*Elastic component  (from Eq 17 in Adhikari et al, GMD 2015): */
		int index=reCast<int,IssmDouble>(alpha/PI*(M-1));
		U_elastic[i] += U_elastic_precomputed[index];
		N_elastic[i] += H_elastic_precomputed[index]*N_azim;
		E_elastic[i] += H_elastic_precomputed[index]*E_azim;

		/*Add all components to the pUp solution vectors:*/
		U_values[i]+=3*rho_ice/rho_earth*area/(4*PI*pow(earth_radius,2))*I*U_elastic[i];
		N_values[i]+=3*rho_ice/rho_earth*area/(4*PI*pow(earth_radius,2))*I*N_elastic[i];
		E_values[i]+=3*rho_ice/rho_earth*area/(4*PI*pow(earth_radius,2))*I*E_elastic[i];
	}

	pUp->SetValues(gsize,indices,U_values,ADD_VAL);
	pNorth->SetValues(gsize,indices,N_values,ADD_VAL);
	pEast->SetValues(gsize,indices,E_values,ADD_VAL);
	
	/*free ressources:*/
	xDelete<int>(indices); 
	xDelete<IssmDouble>(U_values); xDelete<IssmDouble>(N_values); xDelete<IssmDouble>(E_values);
	xDelete<IssmDouble>(U_elastic); xDelete<IssmDouble>(N_elastic); xDelete<IssmDouble>(E_elastic);

	return;
}
/*}}}*/
void    Tria::EsaGeodetic3D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz,IssmDouble eartharea){ /*{{{*/

	/*diverse:*/
	int gsize;
	bool spherical=true;
	IssmDouble llr_list[NUMVERTICES][3];
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble area;
	IssmDouble I;		//ice/water loading 
	IssmDouble late,longe,re;
	IssmDouble lati,longi,ri;
	IssmDouble rho_ice,rho_earth;
	IssmDouble minlong=400;
	IssmDouble maxlong=-20;

	/*precomputed elastic green functions:*/
	IssmDouble* U_elastic_precomputed = NULL;
	IssmDouble* H_elastic_precomputed = NULL;
	int         M;
	
	/*computation of Green functions:*/
	IssmDouble* U_elastic= NULL;
	IssmDouble* N_elastic= NULL;
	IssmDouble* E_elastic= NULL;
	
	/*optimization:*/
	bool store_green_functions=false;

	/*Compute ice thickness change: */
	Input*	deltathickness_input=inputs->GetInput(EsaDeltathicknessEnum); 
	if (!deltathickness_input)_error_("delta thickness input needed to compute elastic adjustment!");
	deltathickness_input->GetInputAverage(&I);
		
	/*early return if we are not on the (ice) loading point: */ 
	if(I==0) return; 

	/*recover material parameters: */
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	rho_earth=matpar->GetMaterialParameter(MaterialsEarthDensityEnum);

	/*how many dofs are we working with here? */
	this->parameters->FindParam(&gsize,MeshNumberofverticesEnum);

	/*compute area of element:*/
	area=GetAreaSpherical();

	/*element centroid (spherical): */
	/* Where is the centroid of this element?:{{{*/
	::GetVerticesCoordinates(&llr_list[0][0],this->vertices,NUMVERTICES,spherical);

	minlong=400; maxlong=-20;
	for (int i=0;i<NUMVERTICES;i++){
		llr_list[i][0]=(90-llr_list[i][0]);
		if(llr_list[i][1]<0)llr_list[i][1]=180+(180+llr_list[i][1]);
		if(llr_list[i][1]>maxlong)maxlong=llr_list[i][1];
		if(llr_list[i][1]<minlong)minlong=llr_list[i][1];
	}
	if(minlong==0 && maxlong>180){
		if (llr_list[0][1]==0)llr_list[0][1]=360;
		if (llr_list[1][1]==0)llr_list[1][1]=360;
		if (llr_list[2][1]==0)llr_list[2][1]=360;
	}

	// correction at the north pole
	if(llr_list[0][0]==0)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==0)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==0)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	//correction at the south pole
	if(llr_list[0][0]==180)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==180)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==180)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	late=(llr_list[0][0]+llr_list[1][0]+llr_list[2][0])/3.0;
	longe=(llr_list[0][1]+llr_list[1][1]+llr_list[2][1])/3.0;

	late=90-late; 
	if(longe>180)longe=(longe-180)-180;

	late=late/180*PI;
	longe=longe/180*PI;
	/*}}}*/

	/*figure out gravity center of our element (Cartesian): */
	IssmDouble x_element, y_element, z_element; 
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x_element=(xyz_list[0][0]+xyz_list[1][0]+xyz_list[2][0])/3.0;
	y_element=(xyz_list[0][1]+xyz_list[1][1]+xyz_list[2][1])/3.0;
	z_element=(xyz_list[0][2]+xyz_list[1][2]+xyz_list[2][2])/3.0;

	/*recover elastic Green's functions for displacement:*/
	DoubleVecParam* U_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(EsaUElasticEnum)); _assert_(U_parameter);
	DoubleVecParam* H_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(EsaHElasticEnum)); _assert_(H_parameter);
	U_parameter->GetParameterValueByPointer(&U_elastic_precomputed,&M);
	H_parameter->GetParameterValueByPointer(&H_elastic_precomputed,&M);

	/*initialize: */
	U_elastic=xNewZeroInit<IssmDouble>(gsize);
	N_elastic=xNewZeroInit<IssmDouble>(gsize);
	E_elastic=xNewZeroInit<IssmDouble>(gsize);

	int* indices=xNew<int>(gsize);
	IssmDouble* U_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* N_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* E_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble alpha;
	IssmDouble delPhi,delLambda;
	IssmDouble dx, dy, dz, x, y, z; 
	IssmDouble N_azim, E_azim;

	for(int i=0;i<gsize;i++){

		indices[i]=i; 

		/*Compute alpha angle between centroid and current vertex: */
		lati=latitude[i]/180*PI; longi=longitude[i]/180*PI;

		delPhi=fabs(lati-late); delLambda=fabs(longi-longe);
		alpha=2.*asin(sqrt(pow(sin(delPhi/2),2.0)+cos(lati)*cos(late)*pow(sin(delLambda/2),2)));

		/*Compute azimuths, both north and east components: */
		x = xx[i]; y = yy[i]; z = zz[i]; 
		if(latitude[i]==90){
			x=1e-12; y=1e-12; 
		}
		if(latitude[i]==-90){
			x=1e-12; y=1e-12; 
		}
		dx = x_element-x; dy = y_element-y; dz = z_element-z; 
		N_azim = (-z*x*dx-z*y*dy+(pow(x,2)+pow(y,2))*dz) /pow((pow(x,2)+pow(y,2))*(pow(x,2)+pow(y,2)+pow(z,2))*(pow(dx,2)+pow(dy,2)+pow(dz,2)),0.5);
		E_azim = (-y*dx+x*dy) /pow((pow(x,2)+pow(y,2))*(pow(dx,2)+pow(dy,2)+pow(dz,2)),0.5);
		
		/*Elastic component  (from Eq 17 in Adhikari et al, GMD 2015): */
		int index=reCast<int,IssmDouble>(alpha/PI*(M-1));
		U_elastic[i] += U_elastic_precomputed[index];
		N_elastic[i] += H_elastic_precomputed[index]*N_azim;
		E_elastic[i] += H_elastic_precomputed[index]*E_azim;

		/*Add all components to the pUp solution vectors:*/
		U_values[i]+=3*rho_ice/rho_earth*area/eartharea*I*U_elastic[i];
		N_values[i]+=3*rho_ice/rho_earth*area/eartharea*I*N_elastic[i];
		E_values[i]+=3*rho_ice/rho_earth*area/eartharea*I*E_elastic[i];
	}
	pUp->SetValues(gsize,indices,U_values,ADD_VAL);
	pNorth->SetValues(gsize,indices,N_values,ADD_VAL);
	pEast->SetValues(gsize,indices,E_values,ADD_VAL);
	
	/*free ressources:*/
	xDelete<int>(indices); 
	xDelete<IssmDouble>(U_values); xDelete<IssmDouble>(N_values); xDelete<IssmDouble>(E_values);
	xDelete<IssmDouble>(U_elastic); xDelete<IssmDouble>(N_elastic); xDelete<IssmDouble>(E_elastic);

	return;
}
/*}}}*/
#endif
#ifdef _HAVE_SEALEVELRISE_
IssmDouble    Tria::OceanArea(void){ /*{{{*/

	if(IsWaterInElement()) return GetAreaSpherical();
	else return 0;

}
/*}}}*/
IssmDouble Tria::OceanAverage(IssmDouble* Sg){ /*{{{*/

	if(IsWaterInElement()){
		
		IssmDouble area;

		/*Compute area of element:*/
		area=GetAreaSpherical();

		/*Average Sg over vertices:*/
		IssmDouble Sg_avg=0; for(int i=0;i<NUMVERTICES;i++) Sg_avg+=Sg[this->vertices[i]->Sid()]/NUMVERTICES;

		/*return: */
		return area*Sg_avg;
	}
	else return 0;

}
/*}}}*/
void	Tria::SealevelriseMomentOfInertia(IssmDouble* dI_list,IssmDouble* Sg_old,IssmDouble eartharea){/*{{{*/

	/*early return if we are not on an ice cap OR ocean:*/
	if(!(this->inputs->Max(MaskIceLevelsetEnum)<0) && !IsWaterInElement()) return; 

	/*Compute area of element:*/
	IssmDouble area; 
	area=GetAreaSpherical();

	/*Compute lat,long,radius of elemental centroid: */
	bool spherical=true;
	IssmDouble llr_list[NUMVERTICES][3];
	IssmDouble late,longe,re;
	/* Where is the centroid of this element?:{{{*/
	::GetVerticesCoordinates(&llr_list[0][0],this->vertices,NUMVERTICES,spherical);

	IssmDouble minlong=400;
	IssmDouble maxlong=-20;
	for (int i=0;i<NUMVERTICES;i++){
		llr_list[i][0]=(90-llr_list[i][0]);
		if(llr_list[i][1]<0)llr_list[i][1]=180+(180+llr_list[i][1]);
		if(llr_list[i][1]>maxlong)maxlong=llr_list[i][1];
		if(llr_list[i][1]<minlong)minlong=llr_list[i][1];
	}
	if(minlong==0 && maxlong>180){
		if (llr_list[0][1]==0)llr_list[0][1]=360;
		if (llr_list[1][1]==0)llr_list[1][1]=360;
		if (llr_list[2][1]==0)llr_list[2][1]=360;
	}

	// correction at the north pole
	if(llr_list[0][0]==0)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==0)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==0)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	//correction at the south pole
	if(llr_list[0][0]==180)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==180)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==180)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	late=(llr_list[0][0]+llr_list[1][0]+llr_list[2][0])/3.0;
	longe=(llr_list[0][1]+llr_list[1][1]+llr_list[2][1])/3.0;

	late=90-late; 
	if(longe>180)longe=(longe-180)-180;

	late=late/180*PI;
	longe=longe/180*PI;
	/*}}}*/
	re=(llr_list[0][2]+llr_list[1][2]+llr_list[2][2])/3.0;

	if(IsWaterInElement()){
		IssmDouble rho_water, S; 
		
		/*recover material parameters: */
		rho_water=matpar->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	
		/*From Sg_old, recover water sea level rise:*/
		S=0; for(int i=0;i<NUMVERTICES;i++) S+=Sg_old[this->vertices[i]->Sid()]/NUMVERTICES;
		
		/* Perturbation terms for moment of inertia (moi_list): 
		 * computed analytically (see Wu & Peltier, eqs 10 & 32) 
		 * also consistent with my GMD formulation!
		 * ALL in geographic coordinates 
		 * */
		dI_list[0] = -4*PI*(rho_water*S*area)*pow(re,4)*(sin(late)*cos(late)*cos(longe))/eartharea; 
		dI_list[1] = -4*PI*(rho_water*S*area)*pow(re,4)*(sin(late)*cos(late)*sin(longe))/eartharea; 
		dI_list[2] = +4*PI*(rho_water*S*area)*pow(re,4)*(1-pow(sin(late),2))/eartharea; 
	}
	else if(this->inputs->Max(MaskIceLevelsetEnum)<0){
		IssmDouble rho_ice, I; 
		
		/*recover material parameters: */
		rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	
		/*Compute ice thickness change: */
		Input*	deltathickness_input=inputs->GetInput(SealevelriseDeltathicknessEnum); 
		if (!deltathickness_input)_error_("delta thickness input needed to compute sea level rise!");
		deltathickness_input->GetInputAverage(&I);
		
		dI_list[0] = -4*PI*(rho_ice*I*area)*pow(re,4)*(sin(late)*cos(late)*cos(longe))/eartharea; 
		dI_list[1] = -4*PI*(rho_ice*I*area)*pow(re,4)*(sin(late)*cos(late)*sin(longe))/eartharea; 
		dI_list[2] = +4*PI*(rho_ice*I*area)*pow(re,4)*(1-pow(sin(late),2))/eartharea; 
	}
	
	return; 
}/*}}}*/
void    Tria::SealevelriseEustatic(Vector<IssmDouble>* pSgi,IssmDouble* peustatic,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble oceanarea,IssmDouble eartharea){ /*{{{*/

	/*diverse:*/
	int gsize;
	bool spherical=true;
	IssmDouble llr_list[NUMVERTICES][3];
	IssmDouble area;
	IssmDouble I;  //change in ice thickness or water level(Farrel and Clarke, Equ. 4)
	IssmDouble rho;
	IssmDouble late,longe,re;
	IssmDouble lati,longi,ri;

	/*elastic green function:*/
	IssmDouble* G_elastic_precomputed=NULL;
	int         M;

	/*ice properties: */
	IssmDouble rho_ice,rho_water,rho_earth;

	/*Initialize eustatic component: do not skip this step :):*/
	IssmDouble eustatic = 0.;

	/*Computational flags:*/
	bool computerigid = true;
	bool computeelastic= true;
	bool scaleoceanarea= false;
	
	/*early return if we are not on an ice cap:*/
	if(!(this->inputs->Max(MaskIceLevelsetEnum)<0)){
		*peustatic=0; //do not forget to assign this pointer, otherwise, global eustatic will be garbage!
		return;
	}

	/*recover material parameters: */
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	rho_water=matpar->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	rho_earth=matpar->GetMaterialParameter(MaterialsEarthDensityEnum);

	/*recover love numbers and computational flags: */
	this->parameters->FindParam(&computerigid,SealevelriseRigidEnum);
	this->parameters->FindParam(&computeelastic,SealevelriseElasticEnum);
	this->parameters->FindParam(&scaleoceanarea,SealevelriseOceanAreaScalingEnum);

	/*recover elastic green function:*/
	if(computeelastic){
		DoubleVecParam* parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelriseGElasticEnum)); 
		_assert_(parameter);
		parameter->GetParameterValueByPointer(&G_elastic_precomputed,&M);
	}

	/*how many dofs are we working with here? */
	this->parameters->FindParam(&gsize,MeshNumberofverticesEnum);

	/* Where is the centroid of this element?:{{{*/
	
	/*retrieve coordinates: */
	::GetVerticesCoordinates(&llr_list[0][0],this->vertices,NUMVERTICES,spherical);
	
	IssmDouble minlong=400;
	IssmDouble maxlong=-20;
	for (int i=0;i<NUMVERTICES;i++){
		llr_list[i][0]=(90-llr_list[i][0]);
		if(llr_list[i][1]<0)llr_list[i][1]=180+(180+llr_list[i][1]);
		if(llr_list[i][1]>maxlong)maxlong=llr_list[i][1];
		if(llr_list[i][1]<minlong)minlong=llr_list[i][1];
	}
	if(minlong==0 && maxlong>180){
		if (llr_list[0][1]==0)llr_list[0][1]=360;
		if (llr_list[1][1]==0)llr_list[1][1]=360;
		if (llr_list[2][1]==0)llr_list[2][1]=360;
	}
	
	// correction at the north pole
	if(llr_list[0][0]==0)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==0)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==0)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;
			
	//correction at the south pole
	if(llr_list[0][0]==180)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==180)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==180)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	late=(llr_list[0][0]+llr_list[1][0]+llr_list[2][0])/3.0;
	longe=(llr_list[0][1]+llr_list[1][1]+llr_list[2][1])/3.0;

	late=90-late; 
	if(longe>180)longe=(longe-180)-180;

	late=late/180*PI;
	longe=longe/180*PI;
	/*}}}*/

	/*Compute area of element:*/
	area=GetAreaSpherical();

	/*Compute ice thickness change: */
	Input*	deltathickness_input=inputs->GetInput(SealevelriseDeltathicknessEnum); 
	if (!deltathickness_input)_error_("delta thickness input needed to compute sea level rise!");
	deltathickness_input->GetInputAverage(&I);

	/*Compute eustatic compoent:*/
	_assert_(oceanarea>0.);
	if(scaleoceanarea) oceanarea=3.619e+14; // use true ocean area, m^2 
	eustatic += rho_ice*area*I/(oceanarea*rho_water); 

	if(computeelastic | computerigid){
		int* indices=xNew<int>(gsize);
		IssmDouble* values=xNew<IssmDouble>(gsize);
		IssmDouble alpha;
		IssmDouble delPhi,delLambda;
		for(int i=0;i<gsize;i++){
			indices[i]=i;

			IssmDouble G_rigid=0;  //do not remove =0!
			IssmDouble G_elastic=0;  //do not remove =0!

			/*Compute alpha angle between centroid and current vertex : */
			lati=latitude[i]/180*PI; longi=longitude[i]/180*PI;

		   delPhi=fabs(lati-late); delLambda=fabs(longi-longe);
			alpha=2.*asin(sqrt(pow(sin(delPhi/2),2.0)+cos(lati)*cos(late)*pow(sin(delLambda/2),2)));

			//Rigid earth gravitational perturbation:
			if(computerigid)G_rigid=1.0/2.0/sin(alpha/2.0);

			//Elastic component  (from Eq 17 in Adhikari et al, GMD 2015)
			if(computeelastic){
				int index=reCast<int,IssmDouble>(alpha/PI*reCast<IssmDouble,int>(M-1));
				G_elastic += G_elastic_precomputed[index];
			}

			/*Add all components to the pSgi or pSgo solution vectors:*/
			values[i]=3*rho_ice/rho_earth*area/eartharea*I*(G_rigid+G_elastic);
		}
		pSgi->SetValues(gsize,indices,values,ADD_VAL);
		
		/*free ressources:*/
		xDelete<IssmDouble>(values);
		xDelete<int>(indices);
	}
	
	/*Assign output pointer:*/
	_assert_(!xIsNan<IssmDouble>(eustatic));
	_assert_(!xIsInf<IssmDouble>(eustatic));
	*peustatic=eustatic;
	return;
}
/*}}}*/
void    Tria::SealevelriseNonEustatic(Vector<IssmDouble>* pSgo,IssmDouble* Sg_old,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble eartharea){ /*{{{*/

	/*diverse:*/
	int gsize;
	bool spherical=true;
	IssmDouble llr_list[NUMVERTICES][3];
	IssmDouble area;
	IssmDouble S;  //change in water water level(Farrel and Clarke, Equ. 4)
	IssmDouble late,longe;
	IssmDouble lati,longi,ri;
	IssmDouble minlong=400;
	IssmDouble maxlong=-20;

	/*precomputed elastic green functions:*/
	IssmDouble* G_elastic_precomputed = NULL;
	int         M;
	
	/*computation of Green functions:*/
	IssmDouble* G_elastic= NULL;
	IssmDouble* G_rigid= NULL;

	/*optimization:*/
	bool store_green_functions=false;

	/*ice properties: */
	IssmDouble rho_ice,rho_water,rho_earth;

	/*Computational flags:*/
	bool computerigid = true;
	bool computeelastic= true;

	/*early return if we are not on the ocean:*/
	if (!IsWaterInElement())return;

	/*recover computational flags: */
	this->parameters->FindParam(&computerigid,SealevelriseRigidEnum);
	this->parameters->FindParam(&computeelastic,SealevelriseElasticEnum);

	/*early return if rigid or elastic not requested:*/
	if(!computerigid && !computeelastic) return;

	/*recover material parameters: */
	rho_water=matpar->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	rho_earth=matpar->GetMaterialParameter(MaterialsEarthDensityEnum);

	/*how many dofs are we working with here? */
	this->parameters->FindParam(&gsize,MeshNumberofverticesEnum);

	/*From Sg_old, recover water sea level rise:*/
	S=0; for(int i=0;i<NUMVERTICES;i++) S+=Sg_old[this->vertices[i]->Sid()]/NUMVERTICES;

	/*Compute area of element:*/
	area=GetAreaSpherical();

	/* Where is the centroid of this element?:{{{*/
	::GetVerticesCoordinates(&llr_list[0][0],this->vertices,NUMVERTICES,spherical);

	minlong=400; maxlong=-20;
	for (int i=0;i<NUMVERTICES;i++){
		llr_list[i][0]=(90-llr_list[i][0]);
		if(llr_list[i][1]<0)llr_list[i][1]=180+(180+llr_list[i][1]);
		if(llr_list[i][1]>maxlong)maxlong=llr_list[i][1];
		if(llr_list[i][1]<minlong)minlong=llr_list[i][1];
	}
	if(minlong==0 && maxlong>180){
		if (llr_list[0][1]==0)llr_list[0][1]=360;
		if (llr_list[1][1]==0)llr_list[1][1]=360;
		if (llr_list[2][1]==0)llr_list[2][1]=360;
	}

	// correction at the north pole
	if(llr_list[0][0]==0)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==0)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==0)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	//correction at the south pole
	if(llr_list[0][0]==180)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==180)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==180)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	late=(llr_list[0][0]+llr_list[1][0]+llr_list[2][0])/3.0;
	longe=(llr_list[0][1]+llr_list[1][1]+llr_list[2][1])/3.0;

	late=90-late; 
	if(longe>180)longe=(longe-180)-180;

	late=late/180*PI;
	longe=longe/180*PI;
	/*}}}*/
	
	if(computeelastic){
	
		/*recover elastic green function:*/
		DoubleVecParam* parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelriseGElasticEnum)); _assert_(parameter);
		parameter->GetParameterValueByPointer(&G_elastic_precomputed,&M);

		/*initialize G_elastic:*/
		G_elastic=xNewZeroInit<IssmDouble>(gsize);
	}
	if(computerigid) G_rigid=xNewZeroInit<IssmDouble>(gsize);

	int* indices=xNew<int>(gsize);
	IssmDouble* values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble alpha;
	IssmDouble delPhi,delLambda;

	for(int i=0;i<gsize;i++){

		indices[i]=i; 

		/*Compute alpha angle between centroid and current vertex : */
		lati=latitude[i]/180*PI; longi=longitude[i]/180*PI;

		delPhi=fabs(lati-late); delLambda=fabs(longi-longe);
		alpha=2.*asin(sqrt(pow(sin(delPhi/2),2.0)+cos(lati)*cos(late)*pow(sin(delLambda/2),2)));

		/*Rigid earth gravitational perturbation: */
		if(computerigid){ 
			G_rigid[i]=1.0/2.0/sin(alpha/2.0); 
			values[i]+=3*rho_water/rho_earth*area/eartharea*S*G_rigid[i];
		}

		/*Elastic component  (from Eq 17 in Adhikari et al, GMD 2015): */
		if(computeelastic){
			int index=reCast<int,IssmDouble>(alpha/PI*(M-1));
			G_elastic[i] += G_elastic_precomputed[index];
			values[i]+=3*rho_water/rho_earth*area/eartharea*S*G_elastic[i];
		}
	}
	
	pSgo->SetValues(gsize,indices,values,ADD_VAL);

	/*free ressources:*/
	xDelete<IssmDouble>(values);
	xDelete<int>(indices);

	/*Free ressources:*/
	if(computeelastic) xDelete<IssmDouble>(G_elastic);
	if(computerigid) xDelete<IssmDouble>(G_rigid);

	return;
}
/*}}}*/
void    Tria::SealevelriseGeodetic(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* Sg,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz,IssmDouble eartharea){ /*{{{*/

	/*diverse:*/
	int gsize;
	bool spherical=true;
	IssmDouble llr_list[NUMVERTICES][3];
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble area;
	IssmDouble I, S;		//change in relative ice thickness and sea level
	IssmDouble late,longe,re;
	IssmDouble lati,longi,ri;
	IssmDouble rho_ice,rho_water,rho_earth;
	IssmDouble minlong=400;
	IssmDouble maxlong=-20;

	/*precomputed elastic green functions:*/
	IssmDouble* U_elastic_precomputed = NULL;
	IssmDouble* H_elastic_precomputed = NULL;
	int         M;
	
	/*computation of Green functions:*/
	IssmDouble* U_elastic= NULL;
	IssmDouble* N_elastic= NULL;
	IssmDouble* E_elastic= NULL;

	/*optimization:*/
	bool store_green_functions=false;

	/*computational flags:*/
	bool computerigid = true;
	bool computeelastic= true;

	/*early return if we are not on the ocean or on an ice cap:*/
	if(!(this->inputs->Max(MaskIceLevelsetEnum)<0) && !IsWaterInElement()) return; 

	/*recover computational flags: */
	this->parameters->FindParam(&computerigid,SealevelriseRigidEnum);
	this->parameters->FindParam(&computeelastic,SealevelriseElasticEnum);
	
	/*early return if elastic not requested:*/
	if(!computeelastic) return;

	/*recover material parameters: */
	rho_ice=matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	rho_water=matpar->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	rho_earth=matpar->GetMaterialParameter(MaterialsEarthDensityEnum);

	/*how many dofs are we working with here? */
	this->parameters->FindParam(&gsize,MeshNumberofverticesEnum);

	/*compute area of element:*/
	area=GetAreaSpherical();

	/*element centroid (spherical): */
	/* Where is the centroid of this element?:{{{*/
	::GetVerticesCoordinates(&llr_list[0][0],this->vertices,NUMVERTICES,spherical);

	minlong=400; maxlong=-20;
	for (int i=0;i<NUMVERTICES;i++){
		llr_list[i][0]=(90-llr_list[i][0]);
		if(llr_list[i][1]<0)llr_list[i][1]=180+(180+llr_list[i][1]);
		if(llr_list[i][1]>maxlong)maxlong=llr_list[i][1];
		if(llr_list[i][1]<minlong)minlong=llr_list[i][1];
	}
	if(minlong==0 && maxlong>180){
		if (llr_list[0][1]==0)llr_list[0][1]=360;
		if (llr_list[1][1]==0)llr_list[1][1]=360;
		if (llr_list[2][1]==0)llr_list[2][1]=360;
	}

	// correction at the north pole
	if(llr_list[0][0]==0)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==0)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==0)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	//correction at the south pole
	if(llr_list[0][0]==180)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==180)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==180)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	late=(llr_list[0][0]+llr_list[1][0]+llr_list[2][0])/3.0;
	longe=(llr_list[0][1]+llr_list[1][1]+llr_list[2][1])/3.0;

	late=90-late; 
	if(longe>180)longe=(longe-180)-180;

	late=late/180*PI;
	longe=longe/180*PI;
	/*}}}*/

	/*figure out gravity center of our element (Cartesian): */
	IssmDouble x_element, y_element, z_element; 
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x_element=(xyz_list[0][0]+xyz_list[1][0]+xyz_list[2][0])/3.0;
	y_element=(xyz_list[0][1]+xyz_list[1][1]+xyz_list[2][1])/3.0;
	z_element=(xyz_list[0][2]+xyz_list[1][2]+xyz_list[2][2])/3.0;

	/*recover elastic Green's functions for displacement:*/
	DoubleVecParam* U_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelriseUElasticEnum)); _assert_(U_parameter);
	DoubleVecParam* H_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelriseHElasticEnum)); _assert_(H_parameter);
	U_parameter->GetParameterValueByPointer(&U_elastic_precomputed,&M);
	H_parameter->GetParameterValueByPointer(&H_elastic_precomputed,&M);

	/*From Sg, recover water sea level rise:*/
	S=0; for(int i=0;i<NUMVERTICES;i++) S+=Sg[this->vertices[i]->Sid()]/NUMVERTICES;
	
	/*Compute ice thickness change: */
	Input*	deltathickness_input=inputs->GetInput(SealevelriseDeltathicknessEnum); 
	if (!deltathickness_input)_error_("delta thickness input needed to compute sea level rise!");
	deltathickness_input->GetInputAverage(&I);
		
	/*initialize: */
	U_elastic=xNewZeroInit<IssmDouble>(gsize);
	N_elastic=xNewZeroInit<IssmDouble>(gsize);
	E_elastic=xNewZeroInit<IssmDouble>(gsize);

	int* indices=xNew<int>(gsize);
	IssmDouble* U_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* N_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* E_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble alpha;
	IssmDouble delPhi,delLambda;
	IssmDouble dx, dy, dz, x, y, z; 
	IssmDouble N_azim, E_azim;

	for(int i=0;i<gsize;i++){

		indices[i]=i; 

		/*Compute alpha angle between centroid and current vertex: */
		lati=latitude[i]/180*PI; longi=longitude[i]/180*PI;

		delPhi=fabs(lati-late); delLambda=fabs(longi-longe);
		alpha=2.*asin(sqrt(pow(sin(delPhi/2),2.0)+cos(lati)*cos(late)*pow(sin(delLambda/2),2)));

		/*Compute azimuths, both north and east components: */
		x = xx[i]; y = yy[i]; z = zz[i]; 
		if(latitude[i]==90){
			x=1e-12; y=1e-12; 
		}
		if(latitude[i]==-90){
			x=1e-12; y=1e-12; 
		}
		dx = x_element-x; dy = y_element-y; dz = z_element-z; 
		N_azim = (-z*x*dx-z*y*dy+(pow(x,2)+pow(y,2))*dz) /pow((pow(x,2)+pow(y,2))*(pow(x,2)+pow(y,2)+pow(z,2))*(pow(dx,2)+pow(dy,2)+pow(dz,2)),0.5);
		E_azim = (-y*dx+x*dy) /pow((pow(x,2)+pow(y,2))*(pow(dx,2)+pow(dy,2)+pow(dz,2)),0.5);
		
		/*Elastic component  (from Eq 17 in Adhikari et al, GMD 2015): */
		int index=reCast<int,IssmDouble>(alpha/PI*(M-1));
		U_elastic[i] += U_elastic_precomputed[index];
		N_elastic[i] += H_elastic_precomputed[index]*N_azim;
		E_elastic[i] += H_elastic_precomputed[index]*E_azim;

		/*Add all components to the pUp solution vectors:*/
		if(this->inputs->Max(MaskIceLevelsetEnum)<0){
			U_values[i]+=3*rho_ice/rho_earth*area/eartharea*I*U_elastic[i];
			N_values[i]+=3*rho_ice/rho_earth*area/eartharea*I*N_elastic[i];
			E_values[i]+=3*rho_ice/rho_earth*area/eartharea*I*E_elastic[i];
		}
		else if(IsWaterInElement()) {
			U_values[i]+=3*rho_water/rho_earth*area/eartharea*S*U_elastic[i];
			N_values[i]+=3*rho_water/rho_earth*area/eartharea*S*N_elastic[i];
			E_values[i]+=3*rho_water/rho_earth*area/eartharea*S*E_elastic[i];
		}
	}
	pUp->SetValues(gsize,indices,U_values,ADD_VAL);
	pNorth->SetValues(gsize,indices,N_values,ADD_VAL);
	pEast->SetValues(gsize,indices,E_values,ADD_VAL);

	/*free ressources:*/
	xDelete<int>(indices); 
	xDelete<IssmDouble>(U_values); xDelete<IssmDouble>(N_values); xDelete<IssmDouble>(E_values);
	xDelete<IssmDouble>(U_elastic); xDelete<IssmDouble>(N_elastic); xDelete<IssmDouble>(E_elastic);

	return;
}
/*}}}*/
#endif

#ifdef _HAVE_DAKOTA_
void       Tria::InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrows, int ncols, int name, int type){/*{{{*/

	int             i,t,row;
	IssmDouble      time;
	TransientInput *transientinput = NULL;
	IssmDouble      values[3];

	/*Check that name is an element input*/
	if(!IsInput(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	switch(type){

		case VertexEnum:
			/*Create transient input: */
			for(t=0;t<ncols;t++){ //ncols is the number of times

				/*create input values: */
				for(i=0;i<3;i++){
					row=this->vertices[i]->Sid();
					values[i]=matrix[ncols*row+t];
				}

				/*time:*/
				time=matrix[(nrows-1)*ncols+t];

				if(t==0) transientinput=new TransientInput(name);
				transientinput->AddTimeInput(new TriaInput(name,values,P1Enum),time);
				transientinput->Configure(parameters);
			}
			this->inputs->AddInput(transientinput);
			break;

		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}

}
/*}}}*/
void       Tria::InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){/*{{{*/

	int i,j;

	/*Check that name is an element input*/
	if(!IsInput(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	switch(type){

		case VertexEnum:

			/*New TriaInput*/
			IssmDouble values[3];

			/*Get values on the 3 vertices*/
			for (i=0;i<3;i++){
				values[i]=vector[this->vertices[i]->Sid()]; //careful, vector of values here is not parallel distributed, but serial distributed (from a serial Dakota core!)
			}

			/*Branch on the specified type of update: */
			switch(name){
				case ThicknessEnum:
					IssmDouble  thickness[3];
					IssmDouble  thickness_init[3];
					IssmDouble  hydrostatic_ratio[3];
					IssmDouble  surface[3];
					IssmDouble  bed[3];

					/*retrieve inputs: */
					GetInputListOnVertices(&thickness_init[0],ThicknessEnum);
					GetInputListOnVertices(&hydrostatic_ratio[0],GeometryHydrostaticRatioEnum);
					GetInputListOnVertices(&bed[0],BaseEnum);
					GetInputListOnVertices(&surface[0],SurfaceEnum);

					/*build new bed and surface: */
					if (this->IsFloating()){
						/*hydrostatic equilibrium: */
						IssmDouble rho_ice,rho_water,di;
						rho_ice   = this->matpar->GetMaterialParameter(MaterialsRhoIceEnum);
						rho_water = this->matpar->GetMaterialParameter(MaterialsRhoSeawaterEnum);
						di        = rho_ice/rho_water;

						/*build new thickness: */
						for (j=0; j<3; j++) {
							/*  for observed/interpolated/hydrostatic thickness, remove scaling from any hydrostatic thickness  */
							if (hydrostatic_ratio[j] >= 0.)
								thickness[j]=values[j]-(values[j]/thickness_init[j]-1.)*hydrostatic_ratio[j]*surface[j]/(1.-di);
							/*  for minimum thickness, don't scale  */
							else
								thickness[j]=thickness_init[j];

							/*  check the computed thickness and update bed*/
							if (thickness[j] < 0.) thickness[j]=1./(1.-di);
							bed[j]=surface[j]-thickness[j];
						}
					}
					else{
						/*build new thickness: */
						for (j=0; j<3; j++) {
							/*  for observed thickness, use scaled value  */
							if (hydrostatic_ratio[j] >= 0.)
								thickness[j]=values[j];
							/*  for minimum thickness, don't scale  */
							else
								thickness[j]=thickness_init[j];
						}

						/*update bed on grounded ice: */
						for(j=0;j<3;j++)bed[j]=surface[j]-thickness[j];
					}

					/*Add new inputs: */
					this->inputs->AddInput(new TriaInput(ThicknessEnum,thickness,P1Enum));
					this->inputs->AddInput(new TriaInput(BaseEnum,bed,P1Enum));
					this->inputs->AddInput(new TriaInput(SurfaceEnum,surface,P1Enum));

					break;
				case MaterialsRheologyBEnum:
					this->inputs->AddInput(new TriaInput(MaterialsRheologyBbarEnum,values,P1Enum));
					break;
				default:
					this->inputs->AddInput(new TriaInput(name,values,P1Enum));
			}
			break;

		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}

}
/*}}}*/
#endif
