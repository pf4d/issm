#include "./DamageEvolutionAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void DamageEvolutionAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	int finiteelement;
	iomodel->FindConstant(&finiteelement,"md.damage.elementinterp");

	/*Fetch parameters: */
	int stabilization;
	iomodel->FindConstant(&stabilization,"md.damage.stabilization");

	/*Do not add constraints in DG,  they are weakly imposed*/
	if(stabilization!=3){
		IoModelToConstraintsx(constraints,iomodel,"md.damage.spcdamage",DamageEvolutionAnalysisEnum,finiteelement);
	}

	/*FCT, constraints are imposed using penalties*/
	if(stabilization==4){
		constraints->ActivatePenaltyMethod(DamageEvolutionAnalysisEnum);
	}
}/*}}}*/
void DamageEvolutionAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Nothing for now*/

}/*}}}*/
void DamageEvolutionAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	int finiteelement;

	iomodel->FindConstant(&finiteelement,"md.damage.elementinterp");
	::CreateNodes(nodes,iomodel,DamageEvolutionAnalysisEnum,finiteelement);
}/*}}}*/
int  DamageEvolutionAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void DamageEvolutionAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int finiteelement;
	bool   ismovingfront;

	iomodel->FindConstant(&finiteelement,"md.damage.elementinterp");
	iomodel->FindConstant(&ismovingfront,"md.transient.ismovingfront");

	/*Update elements: */
	iomodel->FetchData(1,"md.flowequation.element_equation");
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}
	iomodel->DeleteData(1,"md.flowequation.element_equation");

	/*What input do I need to run my damage evolution model?*/
	iomodel->FetchDataToInput(elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vy",VyEnum);
	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchDataToInput(elements,"md.initialization.vz",VzEnum);
	iomodel->FetchDataToInput(elements,"md.damage.D",DamageDEnum);
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.pressure",PressureEnum);

}/*}}}*/
void DamageEvolutionAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*Intermediaries*/
	int         numoutputs;
	char**      requestedoutputs = NULL;

	/*retrieve some parameters: */
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.law",DamageLawEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.stabilization",DamageStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.maxiter",DamageMaxiterEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.max_damage",DamageMaxDamageEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.damage.elementinterp",DamageElementinterpEnum));

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.damage.requested_outputs");
	parameters->AddObject(new IntParam(DamageEvolutionNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(DamageEvolutionRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.damage.requested_outputs");

	/*Retrieve law dependent parameters: */
	int law;
	iomodel->FindConstant(&law,"md.damage.law");
	if (law==0){
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.stress_threshold",DamageStressThresholdEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.kappa",DamageKappaEnum));
	}
	else if (law>0){
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.c1",DamageC1Enum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.c2",DamageC2Enum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.c3",DamageC3Enum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.c4",DamageC4Enum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.stress_threshold",DamageStressThresholdEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.kappa",DamageKappaEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.healing",DamageHealingEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.equiv_stress",DamageEquivStressEnum));
	}

}/*}}}*/

/*Finite Element Analysis*/
void           DamageEvolutionAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           DamageEvolutionAnalysis::CreateDamageFInput(Element* element){/*{{{*/

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	IssmDouble* f   = xNew<IssmDouble>(numnodes);

	/*Calculate damage evolution source term: */
	for (int i=0;i<numnodes;i++){

		/* healing could be handled here */

		/* no source term; damage handled in stress balance */
		f[i]=0.;
	}

	/*Add input*/
	element->AddInput(DamageFEnum,f,element->GetElementType());
	
	/*Clean up and return*/
	xDelete<IssmDouble>(f);
}/*}}}*/
void           DamageEvolutionAnalysis::CreateDamageFInputExp(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble epsf,stress_threshold,eps0;
	IssmDouble damage,B,n,epseff;
	IssmDouble eps_xx,eps_yy,eps_xy,eps1,eps2,epstmp;
	int domaintype;

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	IssmDouble* f   = xNew<IssmDouble>(numnodes);

	/*retrieve parameters:*/
	element->FindParam(&epsf,DamageC1Enum);
	element->FindParam(&stress_threshold,DamageStressThresholdEnum);
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Compute stress tensor: */
	element->ComputeStrainRate();

	/*retrieve what we need: */
	Input* eps_xx_input  = element->GetInput(StrainRatexxEnum);     _assert_(eps_xx_input);
	Input* eps_xy_input  = element->GetInput(StrainRatexyEnum);     _assert_(eps_xy_input);
	Input* eps_yy_input  = element->GetInput(StrainRateyyEnum);     _assert_(eps_yy_input);
	Input*  n_input=element->GetInput(MaterialsRheologyNEnum); _assert_(n_input);
	Input* damage_input = NULL;
	Input* B_input = NULL;
	if(domaintype==Domain2DhorizontalEnum){
		damage_input = element->GetInput(DamageDbarEnum); 	_assert_(damage_input);
		B_input=element->GetInput(MaterialsRheologyBbarEnum); _assert_(B_input);
	}
	else{
		damage_input = element->GetInput(DamageDEnum);   _assert_(damage_input);
		B_input=element->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
	}

	/*Calculate damage evolution source term: */
	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);
		
		eps_xx_input->GetInputValue(&eps_xx,gauss);
		eps_xy_input->GetInputValue(&eps_xy,gauss);
		eps_yy_input->GetInputValue(&eps_yy,gauss);
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		damage_input->GetInputValue(&damage,gauss);
	
		/*Calculate principal effective strain rates*/
		eps1=(eps_xx+eps_yy)/2.+sqrt(pow((eps_xx-eps_yy)/2.,2)+pow(eps_xy,2));
		eps2=(eps_xx+eps_yy)/2.-sqrt(pow((eps_xx-eps_yy)/2.,2)+pow(eps_xy,2));
		if(fabs(eps2)>fabs(eps1)){epstmp=eps2; eps2=eps1; eps1=epstmp;}

		/*Calculate effective strain rate and threshold strain rate*/
		epseff=1./sqrt(2.)*sqrt(eps1*eps1-eps1*eps2+eps2*eps2);
		eps0=pow(stress_threshold/B,n);

		if(epseff>eps0){
			f[i]=1.-pow(eps0/epseff,1./n)*exp(-(epseff-eps0)/(epsf-eps0))-damage;
		}
		else f[i]=0;
	}

	/*Add input*/
	element->AddInput(DamageFEnum,f,element->GetElementType());
	
	/*Clean up and return*/
	xDelete<IssmDouble>(f);
	delete gauss;
}/*}}}*/
void           DamageEvolutionAnalysis::CreateDamageFInputPralong(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble c1,c2,c3,healing,stress_threshold;
	IssmDouble s_xx,s_xy,s_xz,s_yy,s_yz,s_zz,s1,s2,s3,stmp;
	IssmDouble J2s,Chi,Psi,PosPsi,NegPsi;
	IssmDouble damage,tau_xx,tau_xy,tau_xz,tau_yy,tau_yz,tau_zz,stressMaxPrincipal;
	int equivstress,domaintype,dim;

	/*Fetch number of vertices and allocate output*/
	int numnodes = element->GetNumberOfNodes();
	IssmDouble* f   = xNew<IssmDouble>(numnodes);

	/*retrieve parameters:*/
	element->FindParam(&c1,DamageC1Enum);
	element->FindParam(&c2,DamageC2Enum);
	element->FindParam(&c3,DamageC3Enum);
	element->FindParam(&healing,DamageHealingEnum);
	element->FindParam(&stress_threshold,DamageStressThresholdEnum);
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Get problem dimension*/
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 3; break;
		default: _error_("not implemented");
	}
	/*Compute stress tensor and Stress Max Principal: */
	element->ComputeDeviatoricStressTensor();
	if(dim==3){
		/*Only works in 3d because the pressure is defined*/
		element->StressMaxPrincipalCreateInput();
	}
	/*retrieve what we need: */
	Input* tau_xx_input  = element->GetInput(DeviatoricStressxxEnum);     _assert_(tau_xx_input);
	Input* tau_xy_input  = element->GetInput(DeviatoricStressxyEnum);     _assert_(tau_xy_input);
	Input* tau_yy_input  = element->GetInput(DeviatoricStressyyEnum);     _assert_(tau_yy_input);
	Input* tau_xz_input  = NULL;
	Input* tau_yz_input  = NULL;
	Input* tau_zz_input  = NULL;
	Input* stressMaxPrincipal_input = NULL;
	if(dim==3){
		tau_xz_input  = element->GetInput(DeviatoricStressxzEnum);     _assert_(tau_xz_input);
		tau_yz_input  = element->GetInput(DeviatoricStressyzEnum);     _assert_(tau_yz_input);
		tau_zz_input  = element->GetInput(DeviatoricStresszzEnum);     _assert_(tau_zz_input);
		stressMaxPrincipal_input = element->GetInput(StressMaxPrincipalEnum); _assert_(stressMaxPrincipal_input);
	}
	Input* damage_input = NULL;
	if(domaintype==Domain2DhorizontalEnum){
		damage_input = element->GetInput(DamageDbarEnum); 	_assert_(damage_input);
	}
	else{
		damage_input = element->GetInput(DamageDEnum);   _assert_(damage_input);
	}

	/*retrieve the desired type of equivalent stress*/
	element->FindParam(&equivstress,DamageEquivStressEnum);

	/*Calculate damage evolution source term: */
	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);
		
		damage_input->GetInputValue(&damage,gauss);
		tau_xx_input->GetInputValue(&tau_xx,gauss);
		tau_xy_input->GetInputValue(&tau_xy,gauss);
		tau_yy_input->GetInputValue(&tau_yy,gauss);
		if(dim==3){
			tau_xz_input->GetInputValue(&tau_xz,gauss);
			tau_yz_input->GetInputValue(&tau_yz,gauss);
			tau_zz_input->GetInputValue(&tau_zz,gauss);
		}
		/*Calculate effective stress components*/
		s_xx=tau_xx/(1.-damage);
		s_xy=tau_xy/(1.-damage);
		s_yy=tau_yy/(1.-damage);
		if(dim==3){
			s_xz=tau_xz/(1.-damage);
			s_yz=tau_yz/(1.-damage);
			s_zz=tau_zz/(1.-damage);
		}
		/*Calculate principal effective stresses*/
		if(dim==2){
			s1=(s_xx+s_yy)/2.+sqrt(pow((s_xx-s_yy)/2.,2)+pow(s_xy,2));
			s2=(s_xx+s_yy)/2.-sqrt(pow((s_xx-s_yy)/2.,2)+pow(s_xy,2));
			if(fabs(s2)>fabs(s1)){stmp=s2; s2=s1; s1=stmp;}

			if(equivstress==0){ /* von Mises */
				Chi=sqrt(s1*s1-s1*s2+s2*s2);
			}
			else if(equivstress==1){ /* max principal stress */
				Chi=s1;
			}
			Psi=Chi-stress_threshold;
			NegPsi=max(-Chi,0.); /* healing only for compressive stresses */
			PosPsi=max(Psi,0.);
			f[i]= c1*(pow(PosPsi,c2) - healing*pow(NegPsi,c2))*pow((1./(1.-damage)),c3);
		}
		else{
			if(equivstress==1){/* max principal stress */
				stressMaxPrincipal_input->GetInputValue(&stressMaxPrincipal,gauss);
				Chi=stressMaxPrincipal/(1.-damage);
			}
			else if(equivstress==0){/* von Mises */
				Chi=sqrt(((s_xx-s_yy)*(s_xx-s_yy)+(s_yy-s_zz)*(s_yy-s_zz)+(s_zz-s_xx)*(s_zz-s_xx)+6.*(s_xy*s_xy+s_yz*s_yz+s_xz*s_xz))/2.);
			}
			Psi=Chi-stress_threshold;
			NegPsi=max(-Chi,0.); /* healing only for compressive stresses */
			PosPsi=max(Psi,0.);
			f[i]= c1*(pow(PosPsi,c2) - healing*pow(NegPsi,c2))*pow((1./(1.-damage)),c3);
		}
	}
	/*Add input*/
	element->AddInput(DamageFEnum,f,element->GetElementType());
	
	/*Clean up and return*/
	xDelete<IssmDouble>(f);
	delete gauss;
}/*}}}*/
ElementVector* DamageEvolutionAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* DamageEvolutionAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* DamageEvolutionAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;
	/*Intermediaries*/
	int         domaintype,dim;
	int         stabilization;
	IssmDouble  Jdet,dt,D_scalar,h,hx,hy,hz;
	IssmDouble  vel,vx,vy,vz,dvxdx,dvydy,dvzdz,dvx[3],dvy[3],dvz[3];
	IssmDouble *xyz_list  = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 3; break;
		default: _error_("Not implemented yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    B      = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    Bprime = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    D      = xNewZeroInit<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&stabilization,DamageStabilizationEnum);
	Input* vx_input = element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = element->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){
		vz_input=element->GetInput(VzEnum); _assert_(vz_input);
	}

	if(dim==2) h=element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		
		vx_input->GetInputValue(&vx,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);

		if(dim==3){
			vz_input->GetInputValue(&vz,gauss);
			vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
		}

		D_scalar=gauss->weight*Jdet;
		TripleMultiply(basis,1,numnodes,1,
					&D_scalar,1,1,0,
					basis,1,numnodes,0,
					&Ke->values[0],1);

		GetB(B,element,dim,xyz_list,gauss);
		GetBprime(Bprime,element,dim,xyz_list,gauss);

		dvxdx=dvx[0];
		dvydy=dvy[1];
		if(dim==3) dvzdz=dvz[2];
		D_scalar=dt*gauss->weight*Jdet;

		D[0*dim+0]=D_scalar*dvxdx;
		D[1*dim+1]=D_scalar*dvydy;
		if(dim==3) D[2*dim+2]=D_scalar*dvzdz;

		TripleMultiply(B,dim,numnodes,1,
					D,dim,dim,0,
					B,dim,numnodes,0,
					&Ke->values[0],1);

		D[0*dim+0]=D_scalar*vx;
		D[1*dim+1]=D_scalar*vy;
		if(dim==3) D[2*dim+2]=D_scalar*vz;

		TripleMultiply(B,dim,numnodes,1,
					D,dim,dim,0,
					Bprime,dim,numnodes,0,
					&Ke->values[0],1);

		if(stabilization==2){
			if(dim==3){
				vel=sqrt(vx*vx+vy*vy+vz*vz)+1.e-8;
				D[0*dim+0]=h/(2.0*vel)*vx*vx;
				D[1*dim+0]=h/(2.0*vel)*vy*vx;
				D[2*dim+0]=h/(2.0*vel)*vz*vx;
				D[0*dim+1]=h/(2.0*vel)*vx*vy;
				D[1*dim+1]=h/(2.0*vel)*vy*vy;
				D[2*dim+1]=h/(2.0*vel)*vy*vz;
				D[0*dim+2]=h/(2.0*vel)*vx*vz;
				D[1*dim+2]=h/(2.0*vel)*vy*vz;
				D[2*dim+2]=h/(2.0*vel)*vz*vz;
			}
			else{
				/*Streamline upwinding*/
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				D[0*dim+0]=h/(2.0*vel)*vx*vx;
				D[1*dim+0]=h/(2.0*vel)*vy*vx;
				D[0*dim+1]=h/(2.0*vel)*vx*vy;
				D[1*dim+1]=h/(2.0*vel)*vy*vy;
			}
		}
		else if(stabilization==1){
			if(dim==2){
				vx_input->GetInputAverage(&vx);
				vy_input->GetInputAverage(&vy);
				D[0*dim+0]=h/2.0*fabs(vx);
				D[1*dim+1]=h/2.0*fabs(vy);
			}
			else if(dim==3){ 
				element->ElementSizes(&hx,&hy,&hz);
				vel=sqrt(vx*vx + vy*vy + vz*vz)+1.e-14;
				h=sqrt( pow(hx*vx/vel,2) + pow(hy*vy/vel,2) + pow(hz*vz/vel,2));
				D[0*dim+0]=h/(2.*vel)*fabs(vx*vx);  D[0*dim+1]=h/(2.*vel)*fabs(vx*vy); D[0*dim+2]=h/(2.*vel)*fabs(vx*vz);
				D[1*dim+0]=h/(2.*vel)*fabs(vy*vx);  D[1*dim+1]=h/(2.*vel)*fabs(vy*vy); D[1*dim+2]=h/(2.*vel)*fabs(vy*vz);
				D[2*dim+0]=h/(2.*vel)*fabs(vz*vx);  D[2*dim+1]=h/(2.*vel)*fabs(vz*vy); D[2*dim+2]=h/(2.*vel)*fabs(vz*vz);
			}
		}
		if(stabilization==1 || stabilization==2){
			if(dim==2){
				D[0*dim+0]=D_scalar*D[0*dim+0];
				D[1*dim+0]=D_scalar*D[1*dim+0];
				D[0*dim+1]=D_scalar*D[0*dim+1];
				D[1*dim+1]=D_scalar*D[1*dim+1];
			}
			else if(dim==3){
				D[0*dim+0]=D_scalar*D[0*dim+0];
				D[1*dim+0]=D_scalar*D[1*dim+0];
				D[2*dim+0]=D_scalar*D[2*dim+0];
				D[0*dim+1]=D_scalar*D[0*dim+1];
				D[1*dim+1]=D_scalar*D[1*dim+1];
				D[2*dim+1]=D_scalar*D[2*dim+1];
				D[0*dim+2]=D_scalar*D[0*dim+2];
				D[1*dim+2]=D_scalar*D[1*dim+2];
				D[2*dim+2]=D_scalar*D[2*dim+2];
			}
			TripleMultiply(Bprime,dim,numnodes,1,
						D,dim,dim,0,
						Bprime,dim,numnodes,0,
						&Ke->values[0],1);
		}

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(Bprime);
	xDelete<IssmDouble>(D);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* DamageEvolutionAnalysis::CreatePVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int      domaintype,damagelaw;
	IssmDouble  Jdet,dt;
	IssmDouble  f,damage;
	IssmDouble* xyz_list = NULL;
	/*Get element*/
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&damagelaw,DamageLawEnum);
	switch(damagelaw){
		case 0:
			this->CreateDamageFInput(element);
			break;
		case 1:
			this->CreateDamageFInputPralong(element);
			break;
		case 2:
			this->CreateDamageFInputExp(element);
			break;
		default:
			_error_("not implemented yet");
	}

	Input* damaged_input = NULL;
	Input* damagef_input = element->GetInput(DamageFEnum); _assert_(damagef_input);
	if(domaintype==Domain2DhorizontalEnum){
		damaged_input = element->GetInput(DamageDbarEnum); _assert_(damaged_input);
	}
	else{
		damaged_input = element->GetInput(DamageDEnum); _assert_(damaged_input);
	}


	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		damaged_input->GetInputValue(&damage,gauss);
		damagef_input->GetInputValue(&f,gauss);

		for(int i=0;i<numnodes;i++){
			pe->values[i]+=Jdet*gauss->weight*(damage+dt*f)*basis[i];
		}
	}
	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           DamageEvolutionAnalysis::GetB(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3] where Bi is of size 3*NDOF2. 
	 * For node i, Bi can be expressed in the actual coordinate system
	 * by: 
	 *       Bi=[ N ]
	 *          [ N ]
	 * where N is the finiteelement function for node i.
	 *
	 * We assume B_prog has been allocated already, of size: 2x(NDOF1*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions*/
	IssmDouble* basis=xNew<IssmDouble>(numnodes);
	element->NodalFunctions(basis,gauss);

	/*Build B: */
	for(int i=0;i<numnodes;i++){
		for(int j=0;j<dim;j++){
			B[numnodes*j+i] = basis[i];
		}
	}

	/*Clean-up*/
	xDelete<IssmDouble>(basis);
}/*}}}*/
void           DamageEvolutionAnalysis::GetBprime(IssmDouble* Bprime,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B'  matrix. B'=[B1' B2' B3'] where Bi' is of size 3*NDOF2. 
	 * For node i, Bi' can be expressed in the actual coordinate system
	 * by: 
	 *       Bi_prime=[ dN/dx ]
	 *                [ dN/dy ]
	 * where N is the finiteelement function for node i.
	 *
	 * We assume B' has been allocated already, of size: 3x(NDOF2*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* dbasis=xNew<IssmDouble>(dim*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B': */
	for(int i=0;i<numnodes;i++){
		for(int j=0;j<dim;j++){
			Bprime[numnodes*j+i] = dbasis[j*numnodes+i];
		}
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);

}/*}}}*/
void           DamageEvolutionAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,DamageDbarEnum);
}/*}}}*/
void           DamageEvolutionAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           DamageEvolutionAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int domaintype;
	IssmDouble  max_damage;
	int			*doflist = NULL;

	element->FindParam(&domaintype,DomainTypeEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	element->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* newdamage = xNew<IssmDouble>(numnodes);

	/*Get user-supplied max_damage: */
	element->FindParam(&max_damage,DamageMaxDamageEnum);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		newdamage[i]=solution[doflist[i]];
		/*Check solution*/
		if(xIsNan<IssmDouble>(newdamage[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(newdamage[i])) _error_("Inf found in solution vector");
		/*Enforce D < max_damage and D > 0 */
		if(newdamage[i]>max_damage) newdamage[i]=max_damage;
		else if(newdamage[i]<0.)    newdamage[i]=0.;
	}

	/*Get all inputs and parameters*/
	if(domaintype==Domain2DhorizontalEnum){
		element->AddInput(DamageDbarEnum,newdamage,element->GetElementType());
	}
	else{
		element->AddInput(DamageDEnum,newdamage,element->GetElementType());
	}

	/*Free ressources:*/
	xDelete<IssmDouble>(newdamage);
	xDelete<int>(doflist);
}/*}}}*/
void           DamageEvolutionAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
}/*}}}*/

/*Flux Correction Transport*/
ElementMatrix* DamageEvolutionAnalysis::CreateFctKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	IssmDouble Jdet;
	IssmDouble vx,vy;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int dim      = 2;

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    B      = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    Bprime = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    D      = xNewZeroInit<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vxaverage_input=element->GetInput(VxEnum); _assert_(vxaverage_input);
	Input* vyaverage_input=element->GetInput(VyEnum); _assert_(vyaverage_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		GetB(B,element,dim,xyz_list,gauss);
		GetBprime(Bprime,element,dim,xyz_list,gauss);
		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);

		D[0*dim+0] = -gauss->weight*vx*Jdet;
		D[1*dim+1] = -gauss->weight*vy*Jdet;

		TripleMultiply(B,dim,numnodes,1,
					D,dim,dim,0,
					Bprime,dim,numnodes,0,
					&Ke->values[0],1);

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(Bprime);
	xDelete<IssmDouble>(D);
	delete gauss;
	return Ke;
}/*}}}*/
ElementMatrix* DamageEvolutionAnalysis::CreateMassMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	IssmDouble  D,Jdet;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Me     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		D=gauss->weight*Jdet;
		TripleMultiply(basis,1,numnodes,1,
					&D,1,1,0,
					basis,1,numnodes,0,
					&Me->values[0],1);
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return Me;
}/*}}}*/
void           DamageEvolutionAnalysis::FctKMatrix(Matrix<IssmDouble>** pKff,Matrix<IssmDouble>** pKfs,FemModel* femmodel){/*{{{*/

	/*Output*/
	Matrix<IssmDouble>* Kff = NULL;
	Matrix<IssmDouble>* Kfs = NULL;

	/*Initialize Jacobian Matrix*/
	AllocateSystemMatricesx(&Kff,&Kfs,NULL,NULL,femmodel);

	/*Create and assemble matrix*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element*       element = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ElementMatrix* Ke     = this->CreateFctKMatrix(element);
		if(Ke) Ke->AddToGlobal(Kff,Kfs);
		delete Ke;
	}
	Kff->Assemble();
	Kfs->Assemble();

	/*Assign output pointer*/
	*pKff=Kff;
	if(pKfs){
		*pKfs=Kfs;
	}
	else{
		delete Kfs;
	}
}/*}}}*/
void           DamageEvolutionAnalysis::LumpedMassMatrix(Vector<IssmDouble>** pMlff,FemModel* femmodel){/*{{{*/

	/*Intermediaries*/
	int  configuration_type;

	/*Initialize Lumped mass matrix (actually we just save its diagonal)*/
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	int fsize      = femmodel->nodes->NumberOfDofs(configuration_type,FsetEnum);
	int flocalsize = femmodel->nodes->NumberOfDofsLocal(configuration_type,FsetEnum);
	Vector<IssmDouble>* Mlff = new Vector<IssmDouble>(flocalsize,fsize);

	/*Create and assemble matrix*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element*       element = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ElementMatrix* MLe     = this->CreateMassMatrix(element);
		if(MLe){
			MLe->Lump();
			MLe->AddDiagonalToGlobal(Mlff);
		}
		delete MLe;
	}
	Mlff->Assemble();

	/*Assign output pointer*/
	*pMlff=Mlff;
}/*}}}*/
void           DamageEvolutionAnalysis::MassMatrix(Matrix<IssmDouble>** pMff,FemModel* femmodel){/*{{{*/

	/*Initialize Mass matrix*/
	Matrix<IssmDouble> *Mff = NULL;
	AllocateSystemMatricesx(&Mff,NULL,NULL,NULL,femmodel);

	/*Create and assemble matrix*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element*       element = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ElementMatrix* MLe     = this->CreateMassMatrix(element);
		if(MLe){
			MLe->AddToGlobal(Mff);
		}
		delete MLe;
	}
	Mff->Assemble();

	/*Assign output pointer*/
	*pMff=Mff;
}/*}}}*/
