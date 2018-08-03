/*!\file Element.cpp
 * \brief: implementation of the Element object
 */
/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../classes.h"
#include "../../shared/shared.h"
#include "../../modules/SurfaceMassBalancex/SurfaceMassBalancex.h"
/*}}}*/

/*Constructors/destructor/copy*/
Element::Element(){/*{{{*/
	this->id  = -1;
	this->sid = -1;
	this->inputs     = NULL;
	this->nodes      = NULL;
	this->vertices   = NULL;
	this->material   = NULL;
	this->matpar     = NULL;
	this->inputs     = NULL;
	this->parameters = NULL;
	this->element_type_list=NULL;
}/*}}}*/
Element::~Element(){/*{{{*/
	xDelete<int>(element_type_list);
	delete inputs;
}
/*}}}*/

/*Other*/
void       Element::AddInput(Input* input_in){/*{{{*/

	/*Call inputs method*/
	_assert_(this->inputs);
	this->inputs->AddInput(input_in);
}/*}}}*/
void       Element::ComputeLambdaS(){/*{{{*/

	/*Intermediaries*/
	IssmDouble vx,vy,vz,vmag;
	IssmDouble dvx[3],dvy[3],dvz[3],dvmag[3];
	IssmDouble epseff,epsprime;
	int         dim;
	IssmDouble *xyz_list = NULL;

	/*Retrieve all inputs we will be needing: */
	this->GetVerticesCoordinates(&xyz_list);
	parameters->FindParam(&dim,DomainDimensionEnum);
	Input* vx_input=this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input=NULL;
	if(dim==3){vz_input=this->GetInput(VzEnum); _assert_(vz_input);}

	/*Allocate arrays*/
	int numvertices = this->GetNumberOfVertices();
	IssmDouble* lambdas = xNew<IssmDouble>(numvertices);

	/* Start looping on the number of vertices: */
	Gauss* gauss=this->NewGauss();
	for (int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

		/*Get velocity derivatives in all directions*/
		_assert_(dim>1);
		_assert_(vx_input);
		vx_input->GetInputValue(&vx,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		_assert_(vy_input);
		vy_input->GetInputValue(&vy,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		if(dim==3){
			_assert_(vz_input);
			vz_input->GetInputValue(&vz,gauss);
			vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
		}
		else{
			vz = 0.;
			dvz[0] = 0.; dvz[1] = 0.; dvz[2] = 0.;
			dvx[2]= 0.;
			dvy[2]= 0.;
		}
		/*Calculate velocity magnitude and its derivative*/
		vmag = sqrt(vx*vx+vy*vy+vz*vz);
		if(vmag<1e-12){
			vmag=1e-12;
			dvmag[0]=0;
			dvmag[1]=0;
			dvmag[2]=0;
		}
		else{
			dvmag[0]=1./(2*sqrt(vmag))*(2*vx*dvx[0]+2*vy*dvy[0]+2*vz*dvz[0]);
			dvmag[1]=1./(2*sqrt(vmag))*(2*vx*dvx[1]+2*vy*dvy[1]+2*vz*dvz[1]);
			dvmag[2]=1./(2*sqrt(vmag))*(2*vx*dvx[2]+2*vy*dvy[2]+2*vz*dvz[2]);
		}
		EstarStrainrateQuantities(&epseff,&epsprime,vx,vy,vz,vmag,&dvx[0],&dvy[0],&dvz[0],&dvmag[0]);
		lambdas[iv]=EstarLambdaS(epseff,epsprime);
	}

	/*Add Stress tensor components into inputs*/
	this->AddInput(LambdaSEnum,lambdas,P1Enum);
	
	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(lambdas);

}
/*}}}*/
void       Element::ComputeNewDamage(){/*{{{*/

	IssmDouble *xyz_list=NULL;
	IssmDouble  eps_xx,eps_xy,eps_yy,eps_xz,eps_yz,eps_zz,eps_eff;
	IssmDouble  epsmin=1.e-27;
	IssmDouble  eps_0,kappa,sigma_0,B,D,n,envelopeD;
	int         dim,counter=0;
	IssmDouble  k1,k2,threshold=1.e-12;

	/* Retrieve parameters */
	this->GetVerticesCoordinates(&xyz_list);
	this->ComputeStrainRate();
	parameters->FindParam(&dim,DomainDimensionEnum);
	parameters->FindParam(&kappa,DamageKappaEnum);
	parameters->FindParam(&sigma_0,DamageStressThresholdEnum);

	/* Retrieve inputs */
	Input* eps_xx_input=this->GetInput(StrainRatexxEnum); _assert_(eps_xx_input);
	Input* eps_yy_input=this->GetInput(StrainRateyyEnum); _assert_(eps_yy_input);
	Input* eps_xy_input=this->GetInput(StrainRatexyEnum); _assert_(eps_xy_input);
	Input* eps_xz_input=NULL;
	Input* eps_yz_input=NULL;
	Input* eps_zz_input=NULL;
	if(dim==3){
		eps_xz_input=this->GetInput(StrainRatexzEnum); _assert_(eps_xz_input);
		eps_yz_input=this->GetInput(StrainRateyzEnum); _assert_(eps_yz_input);
		eps_zz_input=this->GetInput(StrainRatezzEnum); _assert_(eps_zz_input);
	}

	/* Fetch number of nodes and allocate output*/
   int numnodes = this->GetNumberOfNodes();
   IssmDouble* newD = xNew<IssmDouble>(numnodes);

	/* Retrieve domain-dependent inputs */
	Input* n_input=this->GetInput(MaterialsRheologyNEnum); _assert_(n_input);
   Input* damage_input = NULL;
   Input* B_input = NULL;
	int domaintype;
   parameters->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype==Domain2DhorizontalEnum){
      damage_input = this->GetInput(DamageDbarEnum);  _assert_(damage_input);
      B_input=this->GetInput(MaterialsRheologyBbarEnum); _assert_(B_input);
   }
   else{
      damage_input = this->GetInput(DamageDEnum);   _assert_(damage_input);
      B_input=this->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
   }
	
	/* Start looping on the number of nodes: */
	Gauss* gauss=this->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(this->GetElementType(),i);

		eps_xx_input->GetInputValue(&eps_xx,gauss);
		eps_yy_input->GetInputValue(&eps_yy,gauss);
		eps_xy_input->GetInputValue(&eps_xy,gauss);
		if(dim==3){
			eps_xz_input->GetInputValue(&eps_xz,gauss);
			eps_yz_input->GetInputValue(&eps_yz,gauss);
			eps_zz_input->GetInputValue(&eps_zz,gauss);
		}
		else{eps_xz=0; eps_yz=0; eps_zz=0;}

		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		eps_eff=sqrt(eps_xx*eps_xx+eps_yy*eps_yy+eps_xy*eps_xy+eps_xz*eps_xz+eps_yz*eps_yz+eps_xx*eps_yy+epsmin*epsmin);

		B_input->GetInputValue(&B,gauss);
      n_input->GetInputValue(&n,gauss);
      damage_input->GetInputValue(&D,gauss);

		/* Compute threshold strain rate from threshold stress */
		eps_0=pow(sigma_0/B,n); 

		if(eps_eff>eps_0){
			/* Compute damage on envelope curve for existing level of effective strain rate */
			envelopeD=1.-pow(eps_0/eps_eff,1./n)*exp(-(eps_eff-eps_0)/(eps_0*(kappa-1.)));

			if(envelopeD>D){
				newD[i]=envelopeD;
			}
			else newD[i]=D;
		}
		else newD[i]=D;
	}

	/* Add new damage input to DamageEnum and NewDamageEnum */
	this->AddInput(NewDamageEnum,newD,this->GetElementType());
	if(domaintype==Domain2DhorizontalEnum){
		this->AddInput(DamageDbarEnum,newD,this->GetElementType());
	}
	else{
		this->AddInput(DamageDEnum,newD,this->GetElementType());
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(newD);
	delete gauss;

}/*}}}*/
void       Element::ComputeStrainRate(){/*{{{*/

	int         dim;
	IssmDouble *xyz_list = NULL;
	IssmDouble  epsilon[6];

	/*Retrieve all inputs we will be needing: */
	this->GetVerticesCoordinates(&xyz_list);
	parameters->FindParam(&dim,DomainDimensionEnum);
	Input* vx_input=this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input=NULL;
	if(dim==3){vz_input=this->GetInput(VzEnum); _assert_(vz_input);}

	/*Allocate arrays*/
	int numvertices = this->GetNumberOfVertices();
	IssmDouble* eps_xx = xNew<IssmDouble>(numvertices);
	IssmDouble* eps_yy = xNew<IssmDouble>(numvertices);
	IssmDouble* eps_zz = xNew<IssmDouble>(numvertices);
	IssmDouble* eps_xy = xNew<IssmDouble>(numvertices);
	IssmDouble* eps_xz = xNew<IssmDouble>(numvertices);
	IssmDouble* eps_yz = xNew<IssmDouble>(numvertices);
	IssmDouble* eps_ef = xNew<IssmDouble>(numvertices);

	/* Start looping on the number of vertices: */
	Gauss* gauss=this->NewGauss();
	for (int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

		/*Compute strain rate viscosity and pressure: */
		if(dim==2)
		 this->StrainRateSSA(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
		else
		 this->StrainRateFS(&epsilon[0],xyz_list,gauss,vx_input,vy_input,vz_input);

		if(dim==2){
			 /* epsilon=[exx,eyy,exy];*/
			eps_xx[iv]=epsilon[0]; 
			eps_yy[iv]=epsilon[1];
			eps_xy[iv]=epsilon[2];
			/* eps_eff^2 = 1/2 ( exx^2 + eyy^2 + 2*exy^2 )*/
			eps_ef[iv] = 1./sqrt(2.)*sqrt(epsilon[0]*epsilon[0] + epsilon[1]*epsilon[1] + 2.*epsilon[2]*epsilon[2]);
		}
		else{
			/*epsilon=[exx eyy ezz exy exz eyz]*/
			eps_xx[iv]=epsilon[0]; 
			eps_yy[iv]=epsilon[1];
			eps_zz[iv]=epsilon[2];
			eps_xy[iv]=epsilon[3]; 
			eps_xz[iv]=epsilon[4];
			eps_yz[iv]=epsilon[5];
			/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
			eps_ef[iv] = sqrt(epsilon[0]*epsilon[0] + epsilon[1]*epsilon[1] + epsilon[3]*epsilon[3] +  epsilon[4]*epsilon[4] + epsilon[5]*epsilon[5] + epsilon[0]*epsilon[1]);
		}
	}

	/*Add Stress tensor components into inputs*/
	this->AddInput(StrainRatexxEnum,eps_xx,P1Enum);
	this->AddInput(StrainRatexyEnum,eps_xy,P1Enum);
	this->AddInput(StrainRatexzEnum,eps_xz,P1Enum);
	this->AddInput(StrainRateyyEnum,eps_yy,P1Enum);
	this->AddInput(StrainRateyzEnum,eps_yz,P1Enum);
	this->AddInput(StrainRatezzEnum,eps_zz,P1Enum);
	this->AddInput(StrainRateeffectiveEnum,eps_ef,P1Enum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(eps_xx);
	xDelete<IssmDouble>(eps_yy);
	xDelete<IssmDouble>(eps_zz);
	xDelete<IssmDouble>(eps_xy);
	xDelete<IssmDouble>(eps_xz);
	xDelete<IssmDouble>(eps_yz);
	xDelete<IssmDouble>(eps_ef);

}
/*}}}*/
void       Element::CoordinateSystemTransform(IssmDouble** ptransform,Node** nodes_list,int numnodes,int* cs_array){/*{{{*/

	int         i,counter;
	int         numdofs   = 0;
	IssmDouble  norm;
	IssmDouble *transform = NULL;
	IssmDouble  coord_system[3][3];

	/*Some checks in debugging mode*/
	_assert_(numnodes && nodes_list);

	/*Get total number of dofs*/
	for(i=0;i<numnodes;i++){
		switch(cs_array[i]){
			case PressureEnum: numdofs+=1; break;
			case XYEnum:       numdofs+=2; break;
			case XYZEnum:      numdofs+=3; break;
			default: _error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Allocate and initialize transform matrix*/
	transform=xNew<IssmDouble>(numdofs*numdofs);
	for(i=0;i<numdofs*numdofs;i++) transform[i]=0.0;

	/*Create transform matrix for all nodes (x,y for 2d and x,y,z for 3d). It is a block matrix
	 *for 3 nodes:

	 *     | T1 0  0 |
	 * Q = | 0  T2 0 |
	 *     | 0  0  T3|
	 *
	 * Where T1 is the transform matrix for node 1. It is a simple copy of the coordinate system
	 * associated to this node*/
	counter=0;
	for(i=0;i<numnodes;i++){
		nodes_list[i]->GetCoordinateSystem(&coord_system[0][0]);
		switch(cs_array[i]){
			case PressureEnum:
				/*DO NOT change anything*/
				transform[(numdofs)*(counter) + counter] = 1.;
				counter+=1;
				break;
			case XYEnum:
				/*We remove the z component, we need to renormalize x and y: x=[x1 x2 0] y=[-x2 x1 0]*/
				norm = sqrt( coord_system[0][0]*coord_system[0][0] + coord_system[1][0]*coord_system[1][0]); _assert_(norm>1.e-4);
				transform[(numdofs)*(counter+0) + counter+0] =   coord_system[0][0]/norm;
				transform[(numdofs)*(counter+0) + counter+1] = - coord_system[1][0]/norm;
				transform[(numdofs)*(counter+1) + counter+0] =   coord_system[1][0]/norm;
				transform[(numdofs)*(counter+1) + counter+1] =   coord_system[0][0]/norm;
				counter+=2;
				break;
			case XYZEnum:
				/*The 3 coordinates are changed (x,y,z)*/
				transform[(numdofs)*(counter+0) + counter+0] = coord_system[0][0];
				transform[(numdofs)*(counter+0) + counter+1] = coord_system[0][1];
				transform[(numdofs)*(counter+0) + counter+2] = coord_system[0][2];
				transform[(numdofs)*(counter+1) + counter+0] = coord_system[1][0];
				transform[(numdofs)*(counter+1) + counter+1] = coord_system[1][1];
				transform[(numdofs)*(counter+1) + counter+2] = coord_system[1][2];
				transform[(numdofs)*(counter+2) + counter+0] = coord_system[2][0];
				transform[(numdofs)*(counter+2) + counter+1] = coord_system[2][1];
				transform[(numdofs)*(counter+2) + counter+2] = coord_system[2][2];
				counter+=3;
				break;
			default:
				_error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Assign output pointer*/
	*ptransform=transform;
}
/*}}}*/
void       Element::DeepEcho(void){/*{{{*/

	_printf_(EnumToStringx(this->ObjectEnum())<<" element:\n");
	_printf_("   id : "<<this->id <<"\n");
	_printf_("   sid: "<<this->sid<<"\n");
	if(vertices){
		int numvertices = this->GetNumberOfVertices();
		for(int i=0;i<numvertices;i++) vertices[i]->Echo();
	}
	else _printf_("vertices = NULL\n");

	if(nodes){
		int numnodes = this->GetNumberOfNodes();
		for(int i=0;i<numnodes;i++) nodes[i]->DeepEcho();
	}
	else _printf_("nodes = NULL\n");

	if (material) material->DeepEcho();
	else _printf_("material = NULL\n");

	if (matpar) matpar->DeepEcho();
	else _printf_("matpar = NULL\n");

	_printf_("   parameters\n");
	if (parameters) parameters->DeepEcho();
	else _printf_("parameters = NULL\n");

	_printf_("   inputs\n");
	if (inputs) inputs->DeepEcho();
	else _printf_("inputs=NULL\n");

	return;
}
/*}}}*/
void       Element::DeleteInput(int input_enum){/*{{{*/

	inputs->DeleteInput(input_enum);

}
/*}}}*/
void       Element::DeleteMaterials(void){/*{{{*/
	delete this->material;
}/*}}}*/
void       Element::Delta18oParameterization(void){/*{{{*/

	/*Are we on the base? If not, return*/
	if(!IsOnBase()) return;

	int        numvertices = this->GetNumberOfVertices();

	int        i;
	IssmDouble* monthlytemperatures=xNew<IssmDouble>(12*numvertices);
	IssmDouble* monthlyprec=xNew<IssmDouble>(12*numvertices);
	IssmDouble* TemperaturesPresentday=xNew<IssmDouble>(12*numvertices);
	IssmDouble* TemperaturesLgm=xNew<IssmDouble>(12*numvertices);
	IssmDouble* PrecipitationsPresentday=xNew<IssmDouble>(12*numvertices);
	IssmDouble* tmp=xNew<IssmDouble>(numvertices);
	IssmDouble Delta18oPresent,Delta18oLgm,Delta18oTime;
	IssmDouble Delta18oSurfacePresent,Delta18oSurfaceLgm,Delta18oSurfaceTime;
	IssmDouble time,yts,finaltime,time_yr;

	/*Recover parameters*/
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	this->parameters->FindParam(&finaltime,TimesteppingFinalTimeEnum);
	time_yr=floor(time/yts)*yts;

	/*Recover present day temperature and precipitation*/
	Input* input=this->inputs->GetInput(SmbTemperaturesPresentdayEnum);    _assert_(input);
	Input* input2=this->inputs->GetInput(SmbTemperaturesLgmEnum);          _assert_(input2);
	Input* input3=this->inputs->GetInput(SmbPrecipitationsPresentdayEnum); _assert_(input3);
	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int month=0;month<12;month++){
		for(int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			input->GetInputValue(&TemperaturesPresentday[iv*12+month],gauss,month/12.*yts);
			input2->GetInputValue(&TemperaturesLgm[iv*12+month],gauss,month/12.*yts);
			input3->GetInputValue(&PrecipitationsPresentday[iv*12+month],gauss,month/12.*yts);

			PrecipitationsPresentday[iv*12+month]=PrecipitationsPresentday[iv*12+month]*yts;
		}
	}

	/*Recover delta18o and Delta18oSurface at present day, lgm and at time t*/
	this->parameters->FindParam(&Delta18oPresent,SmbDelta18oEnum,finaltime);
	this->parameters->FindParam(&Delta18oLgm,SmbDelta18oEnum,(finaltime-(21000*yts)));
	this->parameters->FindParam(&Delta18oTime,SmbDelta18oEnum,time);
	this->parameters->FindParam(&Delta18oSurfacePresent,SmbDelta18oSurfaceEnum,finaltime);
	this->parameters->FindParam(&Delta18oSurfaceLgm,SmbDelta18oSurfaceEnum,(finaltime-(21000*yts)));
	this->parameters->FindParam(&Delta18oSurfaceTime,SmbDelta18oSurfaceEnum,time);

	/*Compute the temperature and precipitation*/
	for(int iv=0;iv<numvertices;iv++){
		ComputeDelta18oTemperaturePrecipitation(Delta18oSurfacePresent, Delta18oSurfaceLgm, Delta18oSurfaceTime,
					Delta18oPresent, Delta18oLgm, Delta18oTime,
					&PrecipitationsPresentday[iv*12],
					&TemperaturesLgm[iv*12], &TemperaturesPresentday[iv*12],
					&monthlytemperatures[iv*12], &monthlyprec[iv*12]);
	}

	/*Update inputs*/
	TransientInput* NewTemperatureInput = new TransientInput(SmbMonthlytemperaturesEnum);
	TransientInput* NewPrecipitationInput = new TransientInput(SmbPrecipitationEnum);
	for (int imonth=0;imonth<12;imonth++) {
		for(i=0;i<numvertices;i++) tmp[i]=monthlytemperatures[i*12+imonth];
		switch(this->ObjectEnum()){
			case TriaEnum:  NewTemperatureInput->AddTimeInput(new TriaInput(SmbMonthlytemperaturesEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case PentaEnum: NewTemperatureInput->AddTimeInput(new PentaInput(SmbMonthlytemperaturesEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case TetraEnum: NewTemperatureInput->AddTimeInput(new TetraInput(SmbMonthlytemperaturesEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			default: _error_("Not implemented yet");
		}
		for(i=0;i<numvertices;i++) tmp[i]=monthlyprec[i*12+imonth]/yts;
		switch(this->ObjectEnum()){
			case TriaEnum:  NewPrecipitationInput->AddTimeInput(new TriaInput(SmbPrecipitationEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case PentaEnum: NewPrecipitationInput->AddTimeInput(new PentaInput(SmbPrecipitationEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case TetraEnum: NewPrecipitationInput->AddTimeInput(new TetraInput(SmbPrecipitationEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			default: _error_("Not implemented yet");
		}
	}
	NewTemperatureInput->Configure(this->parameters);
	NewPrecipitationInput->Configure(this->parameters);

	this->inputs->AddInput(NewTemperatureInput);
	this->inputs->AddInput(NewPrecipitationInput);

	switch(this->ObjectEnum()){
		case TriaEnum: break;
		case PentaEnum:
		case TetraEnum:
							this->InputExtrude(SmbMonthlytemperaturesEnum,-1);
							this->InputExtrude(SmbPrecipitationEnum,-1);
							break;
		default: _error_("Not implemented yet");
	}

	/*clean-up*/
	delete gauss;
	xDelete<IssmDouble>(monthlytemperatures);
	xDelete<IssmDouble>(monthlyprec);
	xDelete<IssmDouble>(TemperaturesPresentday);
	xDelete<IssmDouble>(TemperaturesLgm);
	xDelete<IssmDouble>(PrecipitationsPresentday);
	xDelete<IssmDouble>(tmp);

}
/*}}}*/
void       Element::Delta18opdParameterization(void){/*{{{*/
	/*Are we on the base? If not, return*/
	if(!IsOnBase()) return;

	int        numvertices = this->GetNumberOfVertices();

	int        i;
	IssmDouble* monthlytemperatures=xNew<IssmDouble>(12*numvertices);
	IssmDouble* monthlyprec=xNew<IssmDouble>(12*numvertices);
	IssmDouble* TemperaturesPresentday=xNew<IssmDouble>(12*numvertices);
	IssmDouble* PrecipitationsPresentday=xNew<IssmDouble>(12*numvertices);
	IssmDouble* tmp=xNew<IssmDouble>(numvertices);
	IssmDouble Delta18oTime;
	IssmDouble dpermil,f;
	IssmDouble time,yts,time_yr,month;
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	this->parameters->FindParam(&f,SmbFEnum);
	time_yr=floor(time/yts)*yts;

	/*Get some pdd parameters*/
	dpermil=this->matpar->GetMaterialParameter(SmbDpermilEnum);

	/*Recover present day temperature and precipitation*/
	Input*     input=this->inputs->GetInput(SmbTemperaturesPresentdayEnum);    _assert_(input);
	Input*     input2=this->inputs->GetInput(SmbPrecipitationsPresentdayEnum); _assert_(input2);
	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int month=0;month<12;month++) {
		for(int iv=0;iv<numvertices;iv++) {
			gauss->GaussVertex(iv);
			input->GetInputValue(&TemperaturesPresentday[iv*12+month],gauss,month/12.*yts);
			input2->GetInputValue(&PrecipitationsPresentday[iv*12+month],gauss,month/12.*yts);

			PrecipitationsPresentday[iv*12+month]=PrecipitationsPresentday[iv*12+month]*yts;
		}
	}

	/*Recover interpolation parameters at time t*/
	this->parameters->FindParam(&Delta18oTime,SmbDelta18oEnum,time);

	/*Compute the temperature and precipitation*/
	for(int iv=0;iv<numvertices;iv++){
		ComputeD18OTemperaturePrecipitationFromPD(Delta18oTime,dpermil,f,
					&PrecipitationsPresentday[iv*12], &TemperaturesPresentday[iv*12],
					&monthlytemperatures[iv*12], &monthlyprec[iv*12]);
	}

	/*Update inputs*/
	TransientInput* NewTemperatureInput = new TransientInput(SmbMonthlytemperaturesEnum);
	TransientInput* NewPrecipitationInput = new TransientInput(SmbPrecipitationEnum);
	for (int imonth=0;imonth<12;imonth++) {
		for(i=0;i<numvertices;i++) tmp[i]=monthlytemperatures[i*12+imonth];
		switch(this->ObjectEnum()){
			case TriaEnum:  NewTemperatureInput->AddTimeInput(new TriaInput(SmbMonthlytemperaturesEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case PentaEnum: NewTemperatureInput->AddTimeInput(new PentaInput(SmbMonthlytemperaturesEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case TetraEnum: NewTemperatureInput->AddTimeInput(new TetraInput(SmbMonthlytemperaturesEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			default: _error_("Not implemented yet");
		}
		for(i=0;i<numvertices;i++) tmp[i]=monthlyprec[i*12+imonth]/yts;
		switch(this->ObjectEnum()){
			case TriaEnum:  NewPrecipitationInput->AddTimeInput(new TriaInput(SmbPrecipitationEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case PentaEnum: NewPrecipitationInput->AddTimeInput(new PentaInput(SmbPrecipitationEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case TetraEnum: NewPrecipitationInput->AddTimeInput(new TetraInput(SmbPrecipitationEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			default: _error_("Not implemented yet");
		}
	}
	NewTemperatureInput->Configure(this->parameters);
	NewPrecipitationInput->Configure(this->parameters);

	this->inputs->AddInput(NewTemperatureInput);
	this->inputs->AddInput(NewPrecipitationInput);

	switch(this->ObjectEnum()){
		case TriaEnum: break;
		case PentaEnum:
		case TetraEnum:
							this->InputExtrude(SmbMonthlytemperaturesEnum,-1);
							this->InputExtrude(SmbPrecipitationEnum,-1);
							break;
		default: _error_("Not implemented yet");
	}

	/*clean-up*/
	delete gauss;
	xDelete<IssmDouble>(monthlytemperatures);
	xDelete<IssmDouble>(monthlyprec);
	xDelete<IssmDouble>(TemperaturesPresentday);
	xDelete<IssmDouble>(PrecipitationsPresentday);
	xDelete<IssmDouble>(tmp);
	
}
/*}}}*/
IssmDouble Element::Divergence(void){/*{{{*/
	/*Compute element divergence*/

	/*Intermediaries*/
	int        dim;
	IssmDouble Jdet;
	IssmDouble divergence=0.;
	IssmDouble dvx[3],dvy[3],dvz[3];
	IssmDouble *xyz_list = NULL;

	/*Get inputs and parameters*/
	this->FindParam(&dim,DomainDimensionEnum);
	Input* vx_input = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){
		vz_input = this->GetInput(VzEnum); _assert_(vz_input);
	}
	this->GetVerticesCoordinates(&xyz_list);

	Gauss* gauss=this->NewGauss(5);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);
		this->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get strain rate assuming that epsilon has been allocated*/
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		if(dim==2){
			divergence += (dvx[0]+dvy[1])*gauss->weight*Jdet;
		}
		else{
			vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
			divergence += (dvx[0]+dvy[1]+dvz[2])*gauss->weight*Jdet;
		}

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return divergence;
}/*}}}*/
void       Element::dViscositydBFS(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/

	/*Intermediaries*/
	int materialstype;
	IssmDouble dmudB;
	IssmDouble epsilon3d[6];/* epsilon=[exx,eyy,exy,exy,exz,eyz];    */
	IssmDouble epsilon2d[3];/* epsilon=[exx,eyy,exy];    */
	IssmDouble eps_eff;
	IssmDouble eps0=1.e-27;

	if(dim==3){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		this->StrainRateFS(&epsilon3d[0],xyz_list,gauss,vx_input,vy_input,vz_input);
		eps_eff = sqrt(epsilon3d[0]*epsilon3d[0] + epsilon3d[1]*epsilon3d[1] + epsilon3d[3]*epsilon3d[3] +  epsilon3d[4]*epsilon3d[4] + epsilon3d[5]*epsilon3d[5] + epsilon3d[0]*epsilon3d[1]+eps0*eps0);
	}
	else{
		/* eps_eff^2 = 1/2 ( exx^2 + eyy^2 + 2*exy^2 )*/
		this->StrainRateSSA(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = 1./sqrt(2.)*sqrt(epsilon2d[0]*epsilon2d[0] + epsilon2d[1]*epsilon2d[1] + 2.*epsilon2d[2]*epsilon2d[2]);
	}
	/*Get viscosity*/
	materialstype=this->material->ObjectEnum();
	switch(materialstype){
		case MaticeEnum:
			material->GetViscosity_B(&dmudB,eps_eff);
			break;
		case MatestarEnum:
			material->ViscosityBFS(&dmudB,dim,xyz_list,gauss,vx_input,vy_input,vz_input);
			break;
		default: _error_("not supported");
	}

	/*Assign output pointer*/
	*pdmudB=dmudB;

}
/*}}}*/
void       Element::dViscositydBHO(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	int materialstype;
	IssmDouble dmudB;
	IssmDouble epsilon3d[5];/* epsilon=[exx,eyy,exy,exy,exz,eyz];    */
	IssmDouble epsilon2d[2];/* epsilon=[exx,eyy,exy];    */
	IssmDouble eps_eff;
	IssmDouble eps0=1.e-27;

	if(dim==3){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		this->StrainRateHO(&epsilon3d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = sqrt(epsilon3d[0]*epsilon3d[0] + epsilon3d[1]*epsilon3d[1] + epsilon3d[2]*epsilon3d[2] + epsilon3d[3]*epsilon3d[3] +  epsilon3d[4]*epsilon3d[4] + epsilon3d[0]*epsilon3d[1]+eps0*eps0);
	}
	else{
		/* eps_eff^2 = 1/2 ( exx^2 + eyy^2 + 2*exy^2 )*/
		this->StrainRateHO2dvertical(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = 1./sqrt(2.)*sqrt(epsilon2d[0]*epsilon2d[0] + 2.*epsilon2d[1]*epsilon2d[1] + eps0*eps0);
	}
	/*Get viscosity*/
	materialstype=this->material->ObjectEnum();
	switch(materialstype){
		case MaticeEnum:
			material->GetViscosity_B(&dmudB,eps_eff);
			break;
		case MatestarEnum:
			material->ViscosityBHO(&dmudB,dim,xyz_list,gauss,vx_input,vy_input);
			break;
		default: _error_("not supported");
	}

	/*Assign output pointer*/
	*pdmudB=dmudB;

}
/*}}}*/
void       Element::dViscositydBSSA(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	int materialstype;
	IssmDouble dmudB;
	IssmDouble epsilon2d[3];/* epsilon=[exx,eyy,exy];    */
	IssmDouble epsilon1d;   /* epsilon=[exx];    */
	IssmDouble eps_eff;

	 if(dim==2){
		 /* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exx*eyy*/
		 this->StrainRateSSA(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		 eps_eff = sqrt(epsilon2d[0]*epsilon2d[0] + epsilon2d[1]*epsilon2d[1] + epsilon2d[2]*epsilon2d[2] + epsilon2d[0]*epsilon2d[1]);
	 }
	 else{
		 /* eps_eff^2 = 1/2 exx^2*/
		 this->StrainRateSSA1d(&epsilon1d,xyz_list,gauss,vx_input);
		 eps_eff = sqrt(epsilon1d*epsilon1d/2.);
	 }

	/*Get viscosity*/
	materialstype=this->material->ObjectEnum();
	switch(materialstype){
		case MaticeEnum:
			material->GetViscosity_B(&dmudB,eps_eff);
			break;
		case MatestarEnum:
			material->ViscosityBSSA(&dmudB,dim,xyz_list,gauss,vx_input,vy_input);
			break;
		default: _error_("not supported");
	}

	/*Assign output pointer*/
	*pdmudB=dmudB;

}
/*}}}*/
void       Element::dViscositydDSSA(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble dmudB;
	IssmDouble epsilon2d[3];/* epsilon=[exx,eyy,exy];    */
	IssmDouble epsilon1d;   /* epsilon=[exx];    */
	IssmDouble eps_eff;

	if(dim==2){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exx*eyy*/
		this->StrainRateSSA(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = sqrt(epsilon2d[0]*epsilon2d[0] + epsilon2d[1]*epsilon2d[1] + epsilon2d[2]*epsilon2d[2] + epsilon2d[0]*epsilon2d[1]);
	}
	else{
		/* eps_eff^2 = 1/2 exx^2*/
		this->StrainRateSSA1d(&epsilon1d,xyz_list,gauss,vx_input);
		eps_eff = sqrt(epsilon1d*epsilon1d/2.);
	}

	/*Get viscosity*/
	material->GetViscosity_D(&dmudB,eps_eff);

	/*Assign output pointer*/
	*pdmudB=dmudB;

}
/*}}}*/
void       Element::Echo(void){/*{{{*/
	_printf_(EnumToStringx(this->ObjectEnum())<<" element:\n");
	_printf_("   id : "<<this->id <<"\n");
	_printf_("   sid: "<<this->sid<<"\n");
	if(vertices){
		int numvertices = this->GetNumberOfVertices();
		for(int i=0;i<numvertices;i++) vertices[i]->Echo();
	}
	else _printf_("vertices = NULL\n");

	if(nodes){
		int numnodes = this->GetNumberOfNodes();
		for(int i=0;i<numnodes;i++) {
			_printf_("nodes[" << i << "] = " << nodes[i]);	
			nodes[i]->Echo();
		}
	}
	else _printf_("nodes = NULL\n");

	if (material) material->Echo();
	else _printf_("material = NULL\n");

	if (matpar) matpar->Echo();
	else _printf_("matpar = NULL\n");

	_printf_("   parameters\n");
	if (parameters) parameters->Echo();
	else _printf_("parameters = NULL\n");

	_printf_("   inputs\n");
	if (inputs) inputs->Echo();
	else _printf_("inputs=NULL\n");
}
/*}}}*/
IssmDouble Element::EnthalpyDiffusionParameter(IssmDouble enthalpy,IssmDouble pressure){/*{{{*/
	return matpar->GetEnthalpyDiffusionParameter(enthalpy,pressure);
}/*}}}*/
IssmDouble Element::EnthalpyDiffusionParameterVolume(int numvertices,IssmDouble* enthalpy,IssmDouble* pressure){/*{{{*/
	return matpar->GetEnthalpyDiffusionParameterVolume(numvertices,enthalpy,pressure);
}/*}}}*/
void       Element::EnthalpyToThermal(IssmDouble* ptemperature,IssmDouble* pwaterfraction,IssmDouble enthalpy,IssmDouble pressure){/*{{{*/
	matpar->EnthalpyToThermal(ptemperature,pwaterfraction,enthalpy,pressure);
}/*}}}*/
void       Element::FindParam(bool* pvalue,int paramenum){/*{{{*/
	this->parameters->FindParam(pvalue,paramenum);
}/*}}}*/
void       Element::FindParam(int* pvalue,int paramenum){/*{{{*/
	this->parameters->FindParam(pvalue,paramenum);
}/*}}}*/
void       Element::FindParam(IssmDouble* pvalue,int paramenum){/*{{{*/
	this->parameters->FindParam(pvalue,paramenum);
}/*}}}*/
void       Element::FindParam(int** pvalues,int* psize,int paramenum){/*{{{*/
	this->parameters->FindParam(pvalues,psize,paramenum);
}/*}}}*/
void       Element::GetDofList(int** pdoflist,int approximation_enum,int setenum){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = this->GetNumberOfNodes();

	/*First, figure out size of doflist and create it: */
	int numberofdofs=0;
	for(int i=0;i<numnodes;i++) numberofdofs+=nodes[i]->GetNumberOfDofs(approximation_enum,setenum);

	/*Allocate output*/
	int* doflist=xNew<int>(numberofdofs);

	/*Populate: */
	int count=0;
	for(int i=0;i<numnodes;i++){
		nodes[i]->GetDofList(doflist+count,approximation_enum,setenum);
		count+=nodes[i]->GetNumberOfDofs(approximation_enum,setenum);
	}

	/*Assign output pointers:*/
	*pdoflist=doflist;
}
/*}}}*/
void       Element::GetDofListPressure(int** pdoflist,int setenum){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = this->NumberofNodesVelocity();
	int pnumnodes = this->NumberofNodesPressure();

	/*First, figure out size of doflist and create it: */
	int numberofdofs=0;
	for(int i=vnumnodes;i<vnumnodes+pnumnodes;i++) numberofdofs+=nodes[i]->GetNumberOfDofs(FSApproximationEnum,setenum);

	/*Allocate output*/
	int* doflist=xNew<int>(numberofdofs);

	/*Populate: */
	int count=0;
	for(int i=vnumnodes;i<vnumnodes+pnumnodes;i++){
		nodes[i]->GetDofList(doflist+count,FSApproximationEnum,setenum);
		count+=nodes[i]->GetNumberOfDofs(FSApproximationEnum,setenum);
	}

	/*Assign output pointers:*/
	*pdoflist=doflist;
}
/*}}}*/
void       Element::GetDofListVelocity(int** pdoflist,int setenum){/*{{{*/

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = this->NumberofNodesVelocity();

	/*First, figure out size of doflist and create it: */
	int numberofdofs=0;
	for(int i=0;i<numnodes;i++) numberofdofs+=nodes[i]->GetNumberOfDofs(FSvelocityEnum,setenum);

	/*Allocate output*/
	int* doflist=xNew<int>(numberofdofs);

	/*Populate: */
	int count=0;
	for(int i=0;i<numnodes;i++){
		nodes[i]->GetDofList(doflist+count,FSvelocityEnum,setenum);
		count+=nodes[i]->GetNumberOfDofs(FSvelocityEnum,setenum);
	}

	/*Assign output pointers:*/
	*pdoflist=doflist;
}
/*}}}*/
Input*     Element::GetInput(int inputenum){/*{{{*/
	return inputs->GetInput(inputenum);
}/*}}}*/
void       Element::GetInputListOnNodes(IssmDouble* pvalue,int enumtype,IssmDouble defaultvalue){/*{{{*/

	_assert_(pvalue);

	Input *input    = this->GetInput(enumtype);
	int    numnodes = this->GetNumberOfNodes();

	/* Start looping on the number of vertices: */
	if(input){
		Gauss* gauss=this->NewGauss();
		for(int iv=0;iv<numnodes;iv++){
			gauss->GaussNode(this->FiniteElement(),iv);
			input->GetInputValue(&pvalue[iv],gauss);
		}
		delete gauss;
	}
	else{
		for(int iv=0;iv<numnodes;iv++) pvalue[iv]=defaultvalue;
	}
}
/*}}}*/
void       Element::GetInputListOnNodes(IssmDouble* pvalue,int enumtype){/*{{{*/

	_assert_(pvalue);

	int    numnodes = this->GetNumberOfNodes();
	Input *input    = this->GetInput(enumtype);
	if(!input) _error_("Input " << EnumToStringx(enumtype) << " not found in element");

	/* Start looping on the number of vertices: */
	Gauss* gauss=this->NewGauss();
	for(int iv=0;iv<numnodes;iv++){
		gauss->GaussNode(this->FiniteElement(),iv);
		input->GetInputValue(&pvalue[iv],gauss);
	}
	delete gauss;
}
/*}}}*/
void       Element::GetInputListOnNodesVelocity(IssmDouble* pvalue,int enumtype){/*{{{*/

	_assert_(pvalue);

	int    numnodes = this->NumberofNodesVelocity();
	Input *input    = this->GetInput(enumtype);
	if(!input) _error_("Input " << EnumToStringx(enumtype) << " not found in element");

	/* Start looping on the number of vertices: */
	Gauss* gauss=this->NewGauss();
	for(int iv=0;iv<numnodes;iv++){
		gauss->GaussNode(this->VelocityInterpolation(),iv);
		input->GetInputValue(&pvalue[iv],gauss);
	}
	delete gauss;
}
/*}}}*/
void       Element::GetInputListOnVertices(IssmDouble* pvalue,int enumtype){/*{{{*/

	/*Recover input*/
	Input* input=this->GetInput(enumtype);
	if (!input) _error_("Input " << EnumToStringx(enumtype) << " not found in element");

	/*Fetch number vertices for this element*/
	int numvertices = this->GetNumberOfVertices();

	/*Checks in debugging mode*/
	_assert_(pvalue);

	/* Start looping on the number of vertices: */
	Gauss*gauss=this->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);
		input->GetInputValue(&pvalue[iv],gauss);
	}

	/*clean-up*/
	delete gauss;
}
/*}}}*/
void       Element::GetInputListOnVertices(IssmDouble* pvalue,int enumtype,IssmDouble defaultvalue){/*{{{*/

	/*Recover input*/
	Input* input=this->GetInput(enumtype);

	/*Checks in debugging mode*/
	_assert_(pvalue);

	/*Fetch number vertices for this element*/
	int numvertices = this->GetNumberOfVertices();

	/* Start looping on the number of vertices: */
	if (input){
		Gauss* gauss=this->NewGauss();
		for (int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			input->GetInputValue(&pvalue[iv],gauss);
		}
		delete gauss;
	}
	else{
		for(int iv=0;iv<numvertices;iv++) pvalue[iv]=defaultvalue;
	}
}
/*}}}*/
void       Element::GetInputLocalMinMaxOnNodes(IssmDouble* min,IssmDouble* max,IssmDouble* ug){/*{{{*/


	/*Get number of nodes for this element*/
	int numnodes = this->GetNumberOfNodes();

	/*Some checks to avoid segmentation faults*/
	_assert_(ug);
	_assert_(numnodes>0);
	_assert_(nodes);

	/*Get element minimum/maximum*/
	IssmDouble input_min = ug[nodes[0]->GetDof(0,GsetEnum)];
	IssmDouble input_max = input_min;
	for(int i=1;i<numnodes;i++){
		if(ug[nodes[i]->GetDof(0,GsetEnum)] < input_min) input_min = ug[nodes[i]->GetDof(0,GsetEnum)];
		if(ug[nodes[i]->GetDof(0,GsetEnum)] > input_max) input_max = ug[nodes[i]->GetDof(0,GsetEnum)];
	}


	/*Second loop to reassign min and max with local extrema*/
	for(int i=0;i<numnodes;i++){
		if(min[nodes[i]->GetDof(0,GsetEnum)]>input_min) min[nodes[i]->GetDof(0,GsetEnum)] = input_min;
		if(max[nodes[i]->GetDof(0,GsetEnum)]<input_max) max[nodes[i]->GetDof(0,GsetEnum)] = input_max;
	}
}
/*}}}*/
void       Element::GetInputValue(bool* pvalue,int inputenum){/*{{{*/

	Input* input=inputs->GetInput(inputenum);
	if(!input) _error_("Input " << EnumToStringx(inputenum) << " not found in element");
	input->GetInputValue(pvalue);

}/*}}}*/
void       Element::GetInputValue(int* pvalue,int inputenum){/*{{{*/

	Input* input=inputs->GetInput(inputenum);
	if(!input) _error_("Input " << EnumToStringx(inputenum) << " not found in element");
	input->GetInputValue(pvalue);

}/*}}}*/
void       Element::GetInputValue(IssmDouble* pvalue,int inputenum){/*{{{*/

	Input* input=inputs->GetInput(inputenum);
	if(!input) _error_("Input " << EnumToStringx(inputenum) << " not found in element");
	input->GetInputValue(pvalue);

}/*}}}*/
void       Element::GetInputValue(IssmDouble* pvalue,Gauss* gauss,int inputenum){/*{{{*/

	Input* input=inputs->GetInput(inputenum);
	if(!input) _error_("Input " << EnumToStringx(inputenum) << " not found in element");
	input->GetInputValue(pvalue,gauss);

}/*}}}*/
void       Element::GetInputsInterpolations(Vector<IssmDouble>* interpolations){/*{{{*/

	int interpolation;

	/*Go through all inputs and assign interpolation in vector*/
	_assert_(this->inputs);
	for(int i=0;i<this->inputs->Size();i++){
		Input* input=xDynamicCast<Input*>(this->inputs->GetObjectByOffset(i));
		switch(input->ObjectEnum()){
			case BoolInputEnum:
			case DoubleInputEnum:
			case IntInputEnum:
				interpolations->SetValue(input->InstanceEnum(),reCast<IssmDouble>(input->ObjectEnum()),INS_VAL);
				break;
			case TriaInputEnum:
				interpolation = input->GetResultInterpolation();
				interpolations->SetValue(input->InstanceEnum(),interpolation,INS_VAL);
				break;
			default:
				_error_("Input "<<EnumToStringx(input->ObjectEnum())<<" not supported yet");
		}
	}

}/*}}}*/
IssmDouble Element::GetMaterialParameter(int enum_in){/*{{{*/

	_assert_(this->matpar);
	switch(enum_in){ // FIXME: change this to material
		case MaterialsRheologyNEnum:
			return this->material->GetN();
		case MaterialsRheologyBEnum:
			return this->material->GetB();
		case MaterialsRheologyBbarEnum:
			return this->material->GetBbar();
		default:
			return this->matpar->GetMaterialParameter(enum_in);
	}
}/*}}}*/
void       Element::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(nodes);
	int numnodes = this->GetNumberOfNodes();
	for(int i=0;i<numnodes;i++){
		lidlist[i]=nodes[i]->Lid();
	}
}
/*}}}*/
void       Element::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(nodes);
	int numnodes = this->GetNumberOfNodes();
	for(int i=0;i<numnodes;i++){
		sidlist[i]=nodes[i]->Sid();
	}
}
/*}}}*/
void       Element::GetPhi(IssmDouble* phi, IssmDouble*  epsilon, IssmDouble viscosity){/*{{{*/
	/*Compute deformational heating from epsilon and viscosity */

	IssmDouble epsilon_matrix[3][3];
	IssmDouble epsilon_eff;
	IssmDouble epsilon_sqr[3][3];

	/* Build epsilon matrix */
	epsilon_matrix[0][0]=epsilon[0];
	epsilon_matrix[1][0]=epsilon[3];
	epsilon_matrix[2][0]=epsilon[4];
	epsilon_matrix[0][1]=epsilon[3];
	epsilon_matrix[1][1]=epsilon[1];
	epsilon_matrix[2][1]=epsilon[5];
	epsilon_matrix[0][2]=epsilon[4];
	epsilon_matrix[1][2]=epsilon[5];
	epsilon_matrix[2][2]=epsilon[2];

	/* Effective value of epsilon_matrix */
	epsilon_sqr[0][0]=epsilon_matrix[0][0]*epsilon_matrix[0][0];
	epsilon_sqr[1][0]=epsilon_matrix[1][0]*epsilon_matrix[1][0];
	epsilon_sqr[2][0]=epsilon_matrix[2][0]*epsilon_matrix[2][0];
	epsilon_sqr[0][1]=epsilon_matrix[0][1]*epsilon_matrix[0][1];
	epsilon_sqr[1][1]=epsilon_matrix[1][1]*epsilon_matrix[1][1];
	epsilon_sqr[2][1]=epsilon_matrix[2][1]*epsilon_matrix[2][1];
	epsilon_sqr[0][2]=epsilon_matrix[0][2]*epsilon_matrix[0][2];
	epsilon_sqr[1][2]=epsilon_matrix[1][2]*epsilon_matrix[1][2];
	epsilon_sqr[2][2]=epsilon_matrix[2][2]*epsilon_matrix[2][2];
	epsilon_eff=1/sqrt(2.)*sqrt(epsilon_sqr[0][0]+epsilon_sqr[0][1]+ epsilon_sqr[0][2]+ epsilon_sqr[1][0]+ epsilon_sqr[1][1]+ epsilon_sqr[1][2]+ epsilon_sqr[2][0]+ epsilon_sqr[2][1]+ epsilon_sqr[2][2]);

	/*Phi = Tr(sigma * eps) 
	 *    = Tr(sigma'* eps)
	 *    = 2 * eps_eff * sigma'_eff
	 *    = 4 * mu * eps_eff ^2*/
	*phi=4.*epsilon_eff*epsilon_eff*viscosity;
}
/*}}}*/
/* void       Element::GetVectorFromInputs(Vector<IssmDouble>* vector,int input_enum){/\*{{{*\/ */

/* 	/\*Fetch number vertices for this element and allocate arrays*\/ */
/* 	int numvertices = this->GetNumberOfVertices(); */
/* 	int*        vertexpidlist = xNew<int>(numvertices); */
/* 	IssmDouble* values        = xNew<IssmDouble>(numvertices); */

/* 	/\*Fill in values*\/ */
/* 	this->GetVertexPidList(vertexpidlist); */
/* 	this->GetInputListOnVertices(values,input_enum); */
/* 	vector->SetValues(numvertices,vertexpidlist,values,INS_VAL); */

/* 	/\*Clean up*\/ */
/* 	xDelete<int>(vertexpidlist); */
/* 	xDelete<IssmDouble>(values); */

/* } */
/* /\*}}}*\/ */
void       Element::GetVectorFromInputs(Vector<IssmDouble>* vector,int input_enum,int type){/*{{{*/

	/*Fetch number vertices for this element and allocate arrays*/
	int         numvertices = this->GetNumberOfVertices();
	int         numnodes    = this->GetNumberOfNodes();
	int*        doflist     = NULL;
	IssmDouble  value;
	IssmDouble* values      = NULL;
	Input*      input       = NULL;

	switch(type){
		case ElementSIdEnum:
			input=inputs->GetInput(input_enum); _assert_(input);
			input->GetInputAverage(&value);
			vector->SetValue(this->sid,value,INS_VAL);
			break;
		case VertexPIdEnum:
			doflist = xNew<int>(numvertices);
			values = xNew<IssmDouble>(numvertices);
			/*Fill in values*/
			this->GetVertexPidList(doflist);
			this->GetInputListOnVertices(values,input_enum);
			vector->SetValues(numvertices,doflist,values,INS_VAL);
			break;
		case VertexSIdEnum:
			doflist = xNew<int>(numvertices);
			values = xNew<IssmDouble>(numvertices);
			/*Fill in values*/
			this->GetVerticesSidList(doflist);
			this->GetInputListOnVertices(values,input_enum);
			vector->SetValues(numvertices,doflist,values,INS_VAL);
			break;
		case NodesEnum:
			doflist = xNew<int>(numnodes);
			values = xNew<IssmDouble>(numnodes);
			/*Fill in values*/
			this->GetInputListOnNodes(values,input_enum);
			this->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
			vector->SetValues(numnodes,doflist,values,INS_VAL);
			break;
		case NodeSIdEnum:
			doflist = xNew<int>(numnodes);
			values = xNew<IssmDouble>(numnodes);
			/*Fill in values*/
			this->GetNodesSidList(doflist);
			this->GetInputListOnNodes(values,input_enum);
			vector->SetValues(numnodes,doflist,values,INS_VAL);
			break;
		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}

	/*Clean up*/
	xDelete<int>(doflist);
	xDelete<IssmDouble>(values);

}
/*}}}*/
void       Element::GetVertexPidList(int* pidlist){/*{{{*/

	int numvertices = this->GetNumberOfVertices();
	for(int i=0;i<numvertices;i++) pidlist[i]=vertices[i]->Pid();

}
/*}}}*/
void       Element::GetVerticesConnectivityList(int* connectivity){/*{{{*/

	int numvertices = this->GetNumberOfVertices();
	for(int i=0;i<numvertices;i++) connectivity[i]=this->vertices[i]->Connectivity();
}
/*}}}*/
void       Element::GetVerticesCoordinates(IssmDouble** pxyz_list){/*{{{*/

	int         numvertices = this->GetNumberOfVertices();
	IssmDouble* xyz_list    = xNew<IssmDouble>(numvertices*3);
	::GetVerticesCoordinates(xyz_list,this->vertices,numvertices);

	*pxyz_list = xyz_list;

}/*}}}*/
void       Element::GetVerticesSidList(int* sidlist){/*{{{*/

	int numvertices = this->GetNumberOfVertices();
	for(int i=0;i<numvertices;i++) sidlist[i]=this->vertices[i]->Sid();
}
/*}}}*/
IssmDouble Element::GetXcoord(IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	/*output*/
	IssmDouble x;

	/*Create list of x*/
	int         numvertices = this->GetNumberOfVertices();
	IssmDouble* x_list      = xNew<IssmDouble>(numvertices);

	for(int i=0;i<numvertices;i++) x_list[i]=xyz_list[i*3+0];
	ValueP1OnGauss(&x,x_list,gauss);

	xDelete<IssmDouble>(x_list);
	return x;
}/*}}}*/
IssmDouble Element::GetYcoord(IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	/*output*/
	IssmDouble y;

	/*Create list of y*/
	int         numvertices = this->GetNumberOfVertices();
	IssmDouble* y_list      = xNew<IssmDouble>(numvertices);

	for(int i=0;i<numvertices;i++) y_list[i]=xyz_list[i*3+1];
	ValueP1OnGauss(&y,y_list,gauss);

	xDelete<IssmDouble>(y_list);
	return y;
}/*}}}*/
IssmDouble Element::GetZcoord(IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	/*output*/
	IssmDouble z;

	/*Create list of z*/
	int         numvertices = this->GetNumberOfVertices();
	IssmDouble* z_list      = xNew<IssmDouble>(numvertices);

	for(int i=0;i<numvertices;i++) z_list[i]=xyz_list[i*3+2];
	ValueP1OnGauss(&z,z_list,gauss);

	xDelete<IssmDouble>(z_list);
	return z;
}/*}}}*/
void       Element::GradientIndexing(int* indexing,int control_index,bool onsid){/*{{{*/

	/*Get number of controls*/
	int num_controls;
	parameters->FindParam(&num_controls,InversionNumControlParametersEnum);

	/*Get number of vertices*/
	int numvertices = this->GetNumberOfVertices();

	/*get gradient indices*/
	if(onsid){
		for(int i=0;i<numvertices;i++){
			indexing[i]=num_controls*this->vertices[i]->Sid() + control_index;
		}
	}
	else{
		for(int i=0;i<numvertices;i++){
			indexing[i]=num_controls*this->vertices[i]->Pid() + control_index;
		}
	}

}
/*}}}*/
bool       Element::HasNodeOnBase(){/*{{{*/
	return (this->inputs->Max(MeshVertexonbaseEnum)>0.);
}/*}}}*/
bool       Element::HasNodeOnSurface(){/*{{{*/
	return (this->inputs->Max(MeshVertexonsurfaceEnum)>0.);
}/*}}}*/
int        Element::Id(){/*{{{*/

	return this->id;

}
/*}}}*/
void       Element::InputChangeName(int original_enum,int new_enum){/*{{{*/
	this->inputs->ChangeEnum(original_enum,new_enum);
}
/*}}}*/
void       Element::InputCreate(IssmDouble* vector,IoModel* iomodel,int M,int N,int vector_type,int vector_enum,int code){/*{{{*/
    
    /*Intermediaries*/
    int        i,t;
    IssmDouble time;
    
    /*Branch on type of vector: nodal or elementary: */
    if(vector_type==1){ //nodal vector
        
        int         numvertices = this->GetNumberOfVertices();
        int        *vertexids   = xNew<int>(numvertices);
        IssmDouble *values      = xNew<IssmDouble>(numvertices);
        
        /*Recover vertices ids needed to initialize inputs*/
        _assert_(iomodel->elements);
        for(i=0;i<numvertices;i++){
            vertexids[i]=reCast<int>(iomodel->elements[numvertices*this->Sid()+i]); //ids for vertices are in the elements array from Matlab
        }
        
        /*Are we in transient or static? */
        if(M==iomodel->numberofvertices){
            for(i=0;i<numvertices;i++) values[i]=vector[vertexids[i]-1];
            this->AddInput(vector_enum,values,P1Enum);
        }
        else if(M==iomodel->numberofvertices+1){
            /*create transient input: */
            IssmDouble* times = xNew<IssmDouble>(N);
            for(t=0;t<N;t++) times[t] = vector[(M-1)*N+t];
            TransientInput* transientinput=new TransientInput(vector_enum,times,N);
            for(t=0;t<N;t++){
                for(i=0;i<numvertices;i++) values[i]=vector[N*(vertexids[i]-1)+t];
                switch(this->ObjectEnum()){
                    case TriaEnum:  transientinput->AddTimeInput(new TriaInput( vector_enum,values,P1Enum)); break;
                    case PentaEnum: transientinput->AddTimeInput(new PentaInput(vector_enum,values,P1Enum)); break;
                    case TetraEnum: transientinput->AddTimeInput(new TetraInput(vector_enum,values,P1Enum)); break;
                    default: _error_("Not implemented yet");
                }
            }
            this->inputs->AddInput(transientinput);
            xDelete<IssmDouble>(times);
        }
        else _error_("nodal vector is either numberofvertices or numberofvertices+1 long. Field provided (" << EnumToStringx(vector_enum) << ") is " << M << " long");
        
        xDelete<IssmDouble>(values);
        xDelete<int>(vertexids);
    }
    else if(vector_type==2){ //element vector
        
        IssmDouble value;
        
        /*Are we in transient or static? */
        if(M==iomodel->numberofelements){
            if (code==5){ //boolean
                this->inputs->AddInput(new BoolInput(vector_enum,reCast<bool>(vector[this->Sid()])));
            }
            else if (code==6){ //integer
                this->inputs->AddInput(new IntInput(vector_enum,reCast<int>(vector[this->Sid()])));
            }
            else if (code==7){ //IssmDouble
                this->inputs->AddInput(new DoubleInput(vector_enum,vector[this->Sid()]));
            }
            else _error_("could not recognize nature of vector from code " << code);
        }
        else if(M==iomodel->numberofelements+1){
            /*create transient input: */
            IssmDouble* times = xNew<IssmDouble>(N);
            for(t=0;t<N;t++) times[t] = vector[(M-1)*N+t];
            TransientInput* transientinput=new TransientInput(vector_enum,times,N);
            TriaInput* bof=NULL;
            for(t=0;t<N;t++){
                value=vector[N*this->Sid()+t];
                switch(this->ObjectEnum()){
                    case TriaEnum:  transientinput->AddTimeInput(new TriaInput( vector_enum,&value,P0Enum)); break;
                    case PentaEnum: transientinput->AddTimeInput(new PentaInput(vector_enum,&value,P0Enum)); break;
                    case TetraEnum: transientinput->AddTimeInput(new TetraInput(vector_enum,&value,P0Enum)); break;
                    default: _error_("Not implemented yet");
                }
            }
            this->inputs->AddInput(transientinput);
            xDelete<IssmDouble>(times);
        }
        else _error_("element vector is either numberofelements or numberofelements+1 long. Field provided (" << EnumToStringx(vector_enum) << ") is " << M << " long");
    }
    else if(vector_type==3){ //element vector
        
        IssmDouble value;
        
        /*For right now we are static */
        if(M==iomodel->numberofelements){
            /*create transient input: */
            IssmDouble* layers = xNewZeroInit<IssmDouble>(N);;
            for(t=0;t<N;t++) layers[t] = vector[N*this->Sid()+t];
            DoubleArrayInput* arrayinput=new DoubleArrayInput(vector_enum,layers,N);
            this->inputs->AddInput(arrayinput);
            xDelete<IssmDouble>(layers);
        }
        else _error_("element vector is either numberofelements or numberofelements+1 long. Field provided (" << EnumToStringx(vector_enum) << ") is " << M << " long");
    }
    else _error_("Cannot add input for vector type " << vector_type << " (not supported)");
}
/*}}}*/
void       Element::InputDuplicate(int original_enum,int new_enum){/*{{{*/

	if(!IsInput(original_enum)) _error_("Enum "<<EnumToStringx(original_enum)<<" is not in IsInput");

	/*Call inputs method*/
	this->inputs->DuplicateInput(original_enum,new_enum);

}
/*}}}*/
void       Element::InputUpdateFromConstant(int constant, int name){/*{{{*/

	/*Check that name is an element input*/
	if(!IsInput(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	/*update input*/
	this->inputs->AddInput(new IntInput(name,constant));
}
/*}}}*/
void       Element::InputUpdateFromConstant(IssmDouble constant, int name){/*{{{*/

	/*Check that name is an element input*/
	if(!IsInput(name)) return;

	/*update input*/
	this->inputs->AddInput(new DoubleInput(name,constant));
}
/*}}}*/
void       Element::InputUpdateFromConstant(bool constant, int name){/*{{{*/

	/*Check that name is an element input*/
	if(!IsInput(name)) return;

	/*update input*/
	this->inputs->AddInput(new BoolInput(name,constant));
}
/*}}}*/
bool       Element::IsFloating(){/*{{{*/

	bool shelf;
	int  migration_style;
	parameters->FindParam(&migration_style,GroundinglineMigrationEnum);

	if(migration_style==SubelementMigrationEnum || migration_style==SubelementMigration2Enum){ //Floating if all nodes are floating
		if(this->inputs->Max(MaskGroundediceLevelsetEnum) <= 0.) shelf=true;
		else shelf=false;
	}
	else if(migration_style==ContactEnum){
		if(this->inputs->Max(MaskGroundediceLevelsetEnum) > 0.) shelf=false;
		else shelf=true;
	}
	else if(migration_style==NoneEnum || migration_style==AggressiveMigrationEnum || migration_style==SoftMigrationEnum || migration_style==GroundingOnlyEnum){ //Floating if all nodes are floating
		if(this->inputs->Min(MaskGroundediceLevelsetEnum) > 0.) shelf=false;
		else shelf=true;
	}
	else _error_("migration_style not implemented yet");

	return shelf;
}/*}}}*/
bool       Element::IsIceInElement(){/*{{{*/
	return (this->inputs->Min(MaskIceLevelsetEnum)<0.);
}
/*}}}*/
bool       Element::IsInput(int name){/*{{{*/
	if (
				name==ThicknessEnum ||
				name==SurfaceEnum ||
				name==BaseEnum ||
				name==BedEnum ||
				name==BalancethicknessThickeningRateEnum ||
				name==BalancethicknessOmegaEnum ||
				name==SigmaNNEnum ||
				name==SurfaceSlopeXEnum ||
				name==SurfaceSlopeYEnum ||
				name==SmbMassBalanceEnum ||
				name==SmbAccumulationEnum ||
				name==SmbRunoffEnum ||
				name==SmbMeltEnum ||
				name==SmbRefreezeEnum ||
				name==SmbEvaporationEnum ||
				name==SmbIsInitializedEnum ||
				name==BasalforcingsGroundediceMeltingRateEnum ||
				name==BasalforcingsFloatingiceMeltingRateEnum ||
				name==BasalforcingsGeothermalfluxEnum ||
				name==SurfaceAreaEnum||
				name==DamageDEnum ||
				name==DamageDbarEnum ||
				name==PressureEnum ||
				name==VxEnum ||
				name==VyEnum ||
				name==VzEnum ||
				name==VxMeshEnum ||
				name==VyMeshEnum ||
				name==VzMeshEnum ||
				name==InversionVxObsEnum ||
				name==InversionVyObsEnum ||
				name==InversionVzObsEnum ||
				name==TemperatureEnum ||
				name==TemperaturePDDEnum ||
				name==EnthalpyEnum ||
				name==EnthalpyPicardEnum ||
				name==WaterfractionEnum||
				name==WatercolumnEnum || 
				name==FrictionCoefficientEnum ||
				name==FrictionAsEnum ||
				name==FrictionEffectivePressureEnum ||
				name==MaskGroundediceLevelsetEnum ||
				name==MaskIceLevelsetEnum ||
				name==IceMaskNodeActivationEnum ||
				name==LevelsetfunctionSlopeXEnum ||
				name==LevelsetfunctionSlopeYEnum ||
				name==LevelsetfunctionPicardEnum ||
				//name==CalvingCalvingrateEnum ||
				name==GradientEnum ||
				name==OldGradientEnum  ||
				name==ConvergedEnum || 
				name==MaterialsRheologyEEnum ||
				name==MaterialsRheologyEbarEnum ||
				name==MaterialsRheologyBEnum ||
				name==MaterialsRheologyBbarEnum ||
				name==MaterialsRheologyNEnum ||
				name==MaterialsRheologyEcEnum ||
				name==MaterialsRheologyEcbarEnum ||
				name==MaterialsRheologyEsEnum ||
				name==MaterialsRheologyEsbarEnum ||
				name==SealevelEnum || 
				name==SealevelUmotionEnum || 
				name==SealevelNmotionEnum || 
				name==SealevelEmotionEnum || 
				name==SealevelAbsoluteEnum || 
				name==SealevelEustaticEnum || 
				name==SealevelriseDeltathicknessEnum || 
				name==EsaUmotionEnum || 
				name==EsaNmotionEnum || 
				name==EsaEmotionEnum || 
				name==EsaStrainratexxEnum ||
				name==EsaStrainratexyEnum || 
				name==EsaStrainrateyyEnum ||
				name==EsaRotationrateEnum ||
				name==EsaDeltathicknessEnum || 
				name==GiaWEnum || 
				name==GiadWdtEnum ||
				name==SedimentHeadEnum ||
				name==EplHeadEnum ||
				name==SedimentHeadOldEnum ||
				name==EplHeadOldEnum ||
				name==StressIntensityFactorEnum ||
				name==StrainRateparallelEnum ||
				name==StrainRateperpendicularEnum ||
				name==HydrologydcEplThicknessOldEnum ||
				name==HydrologydcEplInitialThicknessEnum ||
				name==HydrologydcEplThicknessEnum ||
				name==HydrologydcMaskEplactiveNodeEnum ||
				name==HydrologyHeadEnum ||
	         name==HydrologyHeadOldEnum ||		
				name==StressbalanceConvergenceNumStepsEnum || 
				name==MeshVertexonbaseEnum || 
				name==FrictionPEnum ||
				name==FrictionQEnum ||
				name==FrictionCoefficientcoulombEnum ||
				name==LoadingforceXEnum ||
				name==LoadingforceYEnum ||
				name==VelEnum ||
				name==VxPicardEnum ||
				name==VyPicardEnum

				) {
					return true;
				}
	else return false;
}
/*}}}*/
bool       Element::IsLandInElement(){/*{{{*/
	return (this->inputs->Max(MaskLandLevelsetEnum)>0.);
}
/*}}}*/
bool       Element::IsWaterInElement(){/*{{{*/
	return (this->inputs->Max(MaskOceanLevelsetEnum)>0.);
}
/*}}}*/
void       Element::LinearFloatingiceMeltingRate(){/*{{{*/

	int numvertices      = this->GetNumberOfVertices();
	IssmDouble  deepwaterel,upperwaterel,deepwatermelt;
	IssmDouble* base     = xNew<IssmDouble>(numvertices);
	IssmDouble* values   = xNew<IssmDouble>(numvertices);

	parameters->FindParam(&deepwaterel,BasalforcingsDeepwaterElevationEnum);
	parameters->FindParam(&upperwaterel,BasalforcingsUpperwaterElevationEnum);
	parameters->FindParam(&deepwatermelt,BasalforcingsDeepwaterMeltingRateEnum);

	this->GetInputListOnVertices(base,BaseEnum);
	for(int i=0;i<numvertices;i++){
		if(base[i]>upperwaterel)      values[i]=0;
		else if (base[i]<deepwaterel) values[i]=deepwatermelt;
		else values[i]=deepwatermelt*(base[i]-upperwaterel)/(deepwaterel-upperwaterel);
	}

	this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,values,P1Enum);
	xDelete<IssmDouble>(base);
	xDelete<IssmDouble>(values);

}/*}}}*/
void       Element::MantlePlumeGeothermalFlux(){/*{{{*/

	int numvertices      = this->GetNumberOfVertices();
	IssmDouble  mantleconductivity,nusselt,dtbg,plumeradius,topplumedepth,bottomplumedepth,plumex,plumey;
	IssmDouble  crustthickness,uppercrustthickness,uppercrustheat,lowercrustheat;
	IssmDouble  crustheat,plumeheat,dt,middleplumedepth,a,e,eprime,A0,lambda,Alambda,dAlambda;
	IssmDouble  x,y,z,c;
	IssmDouble* values   = xNew<IssmDouble>(numvertices);
	IssmDouble *xyz_list = NULL;

	parameters->FindParam(&mantleconductivity,BasalforcingsMantleconductivityEnum);
	parameters->FindParam(&nusselt,BasalforcingsNusseltEnum);
	parameters->FindParam(&dtbg,BasalforcingsDtbgEnum);
	parameters->FindParam(&plumeradius,BasalforcingsPlumeradiusEnum);
	parameters->FindParam(&topplumedepth,BasalforcingsTopplumedepthEnum);
	parameters->FindParam(&bottomplumedepth,BasalforcingsBottomplumedepthEnum);
	parameters->FindParam(&plumex,BasalforcingsPlumexEnum);
	parameters->FindParam(&plumey,BasalforcingsPlumeyEnum);
	parameters->FindParam(&crustthickness,BasalforcingsCrustthicknessEnum);
	parameters->FindParam(&uppercrustthickness,BasalforcingsUppercrustthicknessEnum);
	parameters->FindParam(&uppercrustheat,BasalforcingsUppercrustheatEnum);
	parameters->FindParam(&lowercrustheat,BasalforcingsLowercrustheatEnum);

	this->GetVerticesCoordinates(&xyz_list);
	c=plumeradius;
	a=(bottomplumedepth-topplumedepth)/2.;
	e=pow(a*a-c*c,1./2.)/a;
	A0=(1-pow(e,2.))/pow(e,3.)*(1./2.*log((1+e)/(1-e))-e);
	for(int i=0;i<numvertices;i++){
		y=xyz_list[i*3+0]-plumex;
		z=xyz_list[i*3+1]-plumey;
		x=-(a+topplumedepth+crustthickness);
		lambda=(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2)))))/2;
		dAlambda=(-8*a*pow(c,2)*x*(-2*pow(a,2)+2*pow(c,2)+sqrt(2)*sqrt((a-c)*(a+c))*sqrt(pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2))))))*(pow(a,4)*(pow(y,2)+pow(z,2))+pow(c,4)*(pow(y,2)+pow(z,2))+pow(pow(x,2)+pow(y,2)+pow(z,2),2)*(pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2)))))+pow(c,2)*(pow(x,4)-pow(x,2)*(pow(y,2)+pow(z,2))-(pow(y,2)+pow(z,2))*(2*pow(y,2)+2*pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2))))))+pow(a,2)*(-pow(x,4)+pow(x,2)*(pow(y,2)+pow(z,2))+(pow(y,2)+pow(z,2))*(-2*pow(c,2)+2*(pow(y,2)+pow(z,2))+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2))))))))/(sqrt((a-c)*(a+c))*sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2))))*pow(pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2)))),3.5)*pow(-(sqrt(2)*sqrt((a-c)*(a+c)))+sqrt(pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2))))),2)*(sqrt(2)*sqrt((a-c)*(a+c))+sqrt(pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2)+sqrt(pow(-pow(a,2)-pow(c,2)+pow(x,2)+pow(y,2)+pow(z,2),2)+4*(pow(c,2)*pow(x,2)+pow(a,2)*(-pow(c,2)+pow(y,2)+pow(z,2)))))));
		eprime=pow((a*a-plumeradius*plumeradius)/(a*a+lambda),1./2.);
		Alambda=(1.-e*e)/(e*e*e)*(1./2.*log((1.+eprime)/(1.-eprime))-eprime);
		dt=dtbg-(nusselt-1.)/(1.+A0*(nusselt-1.))*(Alambda*dtbg+x*dtbg*dAlambda);
		plumeheat=mantleconductivity*dt;
		crustheat=uppercrustheat*uppercrustthickness+lowercrustheat*(crustthickness-uppercrustthickness);
		values[i]=crustheat+plumeheat;
	}

	this->AddInput(BasalforcingsGeothermalfluxEnum,values,P1Enum);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(values);

}/*}}}*/
void       Element::MarshallElement(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction,int numanalyses){/*{{{*/
	
	_assert_(this);
	if(marshall_direction==MARSHALLING_BACKWARD){
		inputs=new Inputs();
		nodes = NULL;
	}

	MARSHALLING_ENUM(ElementEnum);
	
	MARSHALLING(id);
	MARSHALLING(sid);
	MARSHALLING(element_type);
	MARSHALLING_DYNAMIC(element_type_list,int,numanalyses);
	inputs->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);

}
/*}}}*/
void       Element::MigrateGroundingLine(IssmDouble* phi_ungrounding){/*{{{*/

	int         numvertices = this->GetNumberOfVertices();
	int        i,migration_style;
	IssmDouble bed_hydro,yts;
	IssmDouble rho_water,rho_ice,density;
	IssmDouble* melting = xNew<IssmDouble>(numvertices);
	IssmDouble* phi     = xNew<IssmDouble>(numvertices);
	IssmDouble* h       = xNew<IssmDouble>(numvertices);
	IssmDouble* s       = xNew<IssmDouble>(numvertices);
	IssmDouble* b       = xNew<IssmDouble>(numvertices);
	IssmDouble* r       = xNew<IssmDouble>(numvertices);
	IssmDouble* sl      = xNew<IssmDouble>(numvertices);

	/*Recover info at the vertices: */
	parameters->FindParam(&migration_style,GroundinglineMigrationEnum);
	parameters->FindParam(&yts,ConstantsYtsEnum);
	GetInputListOnVertices(&h[0],ThicknessEnum);
	GetInputListOnVertices(&s[0],SurfaceEnum);
	GetInputListOnVertices(&b[0],BaseEnum);
	GetInputListOnVertices(&r[0],BedEnum);
	GetInputListOnVertices(&sl[0],SealevelEnum);
	GetInputListOnVertices(&phi[0],MaskGroundediceLevelsetEnum);
	rho_water   = matpar->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	rho_ice     = matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	density     = rho_ice/rho_water;

	/*go through vertices, and update inputs, considering them to be TriaVertex type: */
	for(i=0;i<numvertices;i++){
		/* Contact FS*/
		if(migration_style == ContactEnum && phi_ungrounding[vertices[i]->Pid()]<10){
			phi[i]=phi_ungrounding[vertices[i]->Pid()]; 
			if(phi[i]>=0.) b[i]=r[i];
		}
		else if(migration_style == GroundingOnlyEnum && b[i]<r[i]) b[i]=r[i];
		/*Ice shelf: if bed below bathymetry, impose it at the bathymetry and update surface, elso do nothing */
		else if(phi[i]<=0.){
			if(b[i]<=r[i]){ 
				b[i]        = r[i];
				s[i]        = b[i]+h[i];
			}
		}
		/*Ice sheet: if hydrostatic bed above bathymetry, ice sheet starts to unground, elso do nothing */
		/*Change only if AggressiveMigration or if the ice sheet is in contact with the ocean*/
		else{ // phi>0
			bed_hydro=-density*h[i]+sl[i];
			if (bed_hydro>r[i]){
				/*Unground only if the element is connected to the ice shelf*/
				if(migration_style==AggressiveMigrationEnum || migration_style==SubelementMigrationEnum || migration_style==SubelementMigration2Enum){
					s[i]        = (1-density)*h[i]+sl[i];
					b[i]        = -density*h[i]+sl[i];
				}
				else if(migration_style==SoftMigrationEnum && phi_ungrounding[vertices[i]->Pid()]<0.){
					s[i]        = (1-density)*h[i]+sl[i];
					b[i]        = -density*h[i]+sl[i];
				}
				else{
					if(migration_style!=SoftMigrationEnum && migration_style!=ContactEnum && migration_style!=GroundingOnlyEnum) _error_("Error: migration should be Aggressive, Soft, Subelement, Contact or GroundingOnly");
				}
			}
		}
	}

	/*Recalculate phi*/
	for(i=0;i<numvertices;i++){
		if(migration_style==SoftMigrationEnum){
			bed_hydro=-density*h[i]+sl[i];
			if(phi[i]<0. || bed_hydro<=r[i] || phi_ungrounding[vertices[i]->Pid()]<0.){
				phi[i]=h[i]+(r[i]-sl[i])/density;
			}
		}
		else if(migration_style!=ContactEnum) phi[i]=h[i]+(r[i]-sl[i])/density;
		else{
			/*do nothing*/
		}
	}
	this->AddInput(MaskGroundediceLevelsetEnum,&phi[0],P1Enum);

	/*Update inputs*/
	this->AddInput(SurfaceEnum,&s[0],P1Enum);
	this->AddInput(BaseEnum,&b[0],P1Enum);

	/*Delete*/
	xDelete<IssmDouble>(melting);
	xDelete<IssmDouble>(phi);
	xDelete<IssmDouble>(r);
	xDelete<IssmDouble>(b);
	xDelete<IssmDouble>(s);
	xDelete<IssmDouble>(sl);
	xDelete<IssmDouble>(h);

}
/*}}}*/
void       Element::MismipFloatingiceMeltingRate(){/*{{{*/

	int numvertices      = this->GetNumberOfVertices();
	IssmDouble  meltratefactor,thresholdthickness,upperdepthmelt;
	IssmDouble* base     = xNew<IssmDouble>(numvertices);
	IssmDouble* bed      = xNew<IssmDouble>(numvertices);
	IssmDouble* values   = xNew<IssmDouble>(numvertices);

	parameters->FindParam(&meltratefactor,BasalforcingsMeltrateFactorEnum);
	parameters->FindParam(&thresholdthickness,BasalforcingsThresholdThicknessEnum);
	parameters->FindParam(&upperdepthmelt,BasalforcingsUpperdepthMeltEnum);

	this->GetInputListOnVertices(base,BaseEnum);
	this->GetInputListOnVertices(bed,BedEnum);
	for(int i=0;i<numvertices;i++){
		if(base[i]>upperdepthmelt){
			values[i]=0;
		}
		else{
			values[i]=meltratefactor*tanh((base[i]-bed[i])/thresholdthickness)*(upperdepthmelt-base[i]);
		}
	}

	this->AddInput(BasalforcingsFloatingiceMeltingRateEnum,values,P1Enum);
	xDelete<IssmDouble>(base);
	xDelete<IssmDouble>(bed);
	xDelete<IssmDouble>(values);

}/*}}}*/
void       Element::MungsmtpParameterization(void){/*{{{*/
	/*Are we on the base? If not, return*/
	if(!IsOnBase()) return;

	int        numvertices = this->GetNumberOfVertices();

	int        i;
	IssmDouble* monthlytemperatures=xNew<IssmDouble>(12*numvertices);
	IssmDouble* monthlyprec=xNew<IssmDouble>(12*numvertices);
	IssmDouble* TemperaturesPresentday=xNew<IssmDouble>(12*numvertices);
	IssmDouble* TemperaturesLgm=xNew<IssmDouble>(12*numvertices);
	IssmDouble* PrecipitationsPresentday=xNew<IssmDouble>(12*numvertices);
	IssmDouble* PrecipitationsLgm=xNew<IssmDouble>(12*numvertices);
	IssmDouble* tmp=xNew<IssmDouble>(numvertices);
	IssmDouble TdiffTime,PfacTime;
	IssmDouble time,yts,time_yr;
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	time_yr=floor(time/yts)*yts;

	/*Recover present day temperature and precipitation*/
	Input*     input=this->inputs->GetInput(SmbTemperaturesPresentdayEnum);    _assert_(input);
	Input*     input2=this->inputs->GetInput(SmbTemperaturesLgmEnum);          _assert_(input2);
	Input*     input3=this->inputs->GetInput(SmbPrecipitationsPresentdayEnum); _assert_(input3);
	Input*     input4=this->inputs->GetInput(SmbPrecipitationsLgmEnum);        _assert_(input4);
	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int month=0;month<12;month++) {
		for(int iv=0;iv<numvertices;iv++) {
			gauss->GaussVertex(iv);
			input->GetInputValue(&TemperaturesPresentday[iv*12+month],gauss,month/12.*yts);
			input2->GetInputValue(&TemperaturesLgm[iv*12+month],gauss,month/12.*yts);
			input3->GetInputValue(&PrecipitationsPresentday[iv*12+month],gauss,month/12.*yts);
			input4->GetInputValue(&PrecipitationsLgm[iv*12+month],gauss,month/12.*yts);

			PrecipitationsPresentday[iv*12+month]=PrecipitationsPresentday[iv*12+month]*yts;
			PrecipitationsLgm[iv*12+month]=PrecipitationsLgm[iv*12+month]*yts;
		}
	}

	/*Recover interpolation parameters at time t*/
	this->parameters->FindParam(&TdiffTime,SmbTdiffEnum,time);
	this->parameters->FindParam(&PfacTime,SmbPfacEnum,time);

	/*Compute the temperature and precipitation*/
	for(int iv=0;iv<numvertices;iv++){
		ComputeMungsmTemperaturePrecipitation(TdiffTime,PfacTime,
					&PrecipitationsLgm[iv*12],&PrecipitationsPresentday[iv*12],
					&TemperaturesLgm[iv*12], &TemperaturesPresentday[iv*12],
					&monthlytemperatures[iv*12], &monthlyprec[iv*12]);
	}

	/*Update inputs*/
	TransientInput* NewTemperatureInput = new TransientInput(SmbMonthlytemperaturesEnum);
	TransientInput* NewPrecipitationInput = new TransientInput(SmbPrecipitationEnum);
	for (int imonth=0;imonth<12;imonth++) {
		for(i=0;i<numvertices;i++) tmp[i]=monthlytemperatures[i*12+imonth];
		switch(this->ObjectEnum()){
			case TriaEnum:  NewTemperatureInput->AddTimeInput(new TriaInput(SmbMonthlytemperaturesEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case PentaEnum: NewTemperatureInput->AddTimeInput(new PentaInput(SmbMonthlytemperaturesEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case TetraEnum: NewTemperatureInput->AddTimeInput(new TetraInput(SmbMonthlytemperaturesEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			default: _error_("Not implemented yet");
		}
		for(i=0;i<numvertices;i++) tmp[i]=monthlyprec[i*12+imonth]/yts;
		switch(this->ObjectEnum()){
			case TriaEnum:  NewPrecipitationInput->AddTimeInput(new TriaInput(SmbPrecipitationEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case PentaEnum: NewPrecipitationInput->AddTimeInput(new PentaInput(SmbPrecipitationEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			case TetraEnum: NewPrecipitationInput->AddTimeInput(new TetraInput(SmbPrecipitationEnum,&tmp[0],P1Enum),time_yr+imonth/12.*yts); break;
			default: _error_("Not implemented yet");
		}
	}
	NewTemperatureInput->Configure(this->parameters);
	NewPrecipitationInput->Configure(this->parameters);

	this->inputs->AddInput(NewTemperatureInput);
	this->inputs->AddInput(NewPrecipitationInput);

	switch(this->ObjectEnum()){
		case TriaEnum: break;
		case PentaEnum:
		case TetraEnum:
							this->InputExtrude(SmbMonthlytemperaturesEnum,-1);
							this->InputExtrude(SmbPrecipitationEnum,-1);
							break;
		default: _error_("Not implemented yet");
	}

	/*clean-up*/
	delete gauss;
	xDelete<IssmDouble>(monthlytemperatures);
	xDelete<IssmDouble>(monthlyprec);
	xDelete<IssmDouble>(TemperaturesPresentday);
	xDelete<IssmDouble>(TemperaturesLgm);
	xDelete<IssmDouble>(PrecipitationsPresentday);
	xDelete<IssmDouble>(PrecipitationsLgm);
	xDelete<IssmDouble>(tmp);

}
/*}}}*/
ElementMatrix* Element::NewElementMatrix(int approximation_enum){/*{{{*/
	return new ElementMatrix(nodes,this->GetNumberOfNodes(),this->parameters,approximation_enum);
}
/*}}}*/
ElementMatrix* Element::NewElementMatrixCoupling(int number_nodes,int approximation_enum){/*{{{*/
	return new ElementMatrix(nodes,number_nodes,this->parameters,approximation_enum);
}
/*}}}*/
ElementVector* Element::NewElementVector(int approximation_enum){/*{{{*/
	return new ElementVector(nodes,this->GetNumberOfNodes(),this->parameters,approximation_enum);
}
/*}}}*/
void       Element::PositiveDegreeDay(IssmDouble* pdds,IssmDouble* pds,IssmDouble signorm,bool ismungsm){/*{{{*/

	int  numvertices = this->GetNumberOfVertices();

	int        i;
	IssmDouble* agd=xNew<IssmDouble>(numvertices); // surface mass balance
	IssmDouble* melt=xNew<IssmDouble>(numvertices); // surface mass balance
	IssmDouble* accu=xNew<IssmDouble>(numvertices); // surface mass balance
	IssmDouble* monthlytemperatures=xNew<IssmDouble>(12*numvertices);
	IssmDouble* monthlyprec=xNew<IssmDouble>(12*numvertices);
	IssmDouble* yearlytemperatures=xNew<IssmDouble>(numvertices); memset(yearlytemperatures, 0., numvertices*sizeof(IssmDouble));
	IssmDouble* tmp=xNew<IssmDouble>(numvertices);
	IssmDouble* h=xNew<IssmDouble>(numvertices);
	IssmDouble* s=xNew<IssmDouble>(numvertices);
	IssmDouble* s0p=xNew<IssmDouble>(numvertices);
	IssmDouble* s0t=xNew<IssmDouble>(numvertices);
	IssmDouble rho_water,rho_ice,desfac,rlaps,rlapslgm;
	IssmDouble PfacTime,TdiffTime,sealevTime;
	IssmDouble mavg=1./12.; //factor for monthly average

	/*Get material parameters :*/
	rho_water=this->matpar->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	rho_ice=this->matpar->GetMaterialParameter(MaterialsRhoIceEnum);

	/*Get some pdd parameters*/
	desfac=this->matpar->GetMaterialParameter(SmbDesfacEnum);
	rlaps=this->matpar->GetMaterialParameter(SmbRlapsEnum);
	rlapslgm=this->matpar->GetMaterialParameter(SmbRlapslgmEnum);

	/*Recover monthly temperatures and precipitation and compute the yearly mean temperatures*/
	Input*     input=this->inputs->GetInput(SmbMonthlytemperaturesEnum); _assert_(input);
	Input*     input2=this->inputs->GetInput(SmbPrecipitationEnum); _assert_(input2);
	IssmDouble time,yts,time_yr;
	this->parameters->FindParam(&time,TimeEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);
	time_yr=floor(time/yts)*yts;

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for(int month=0;month<12;month++) {
		for(int iv=0;iv<numvertices;iv++) {
			gauss->GaussVertex(iv);
			input->GetInputValue(&monthlytemperatures[iv*12+month],gauss,time_yr+month/12.*yts);
			// yearlytemperatures[iv]=yearlytemperatures[iv]+monthlytemperatures[iv*12+month]*mavg; // Has to be in Kelvin
			monthlytemperatures[iv*12+month]=monthlytemperatures[iv*12+month]-273.15; // conversion from Kelvin to celcius for PDD module
			input2->GetInputValue(&monthlyprec[iv*12+month],gauss,time_yr+month/12.*yts);
			monthlyprec[iv*12+month]=monthlyprec[iv*12+month]*yts;
		}
	}

	/*Recover Pfac, Tdiff and sealev at time t:
	 *     This parameters are used to interpolate the temperature
	 *         and precipitaton between PD and LGM when ismungsm==1 */
	if (ismungsm==1){
		this->parameters->FindParam(&TdiffTime,SmbTdiffEnum,time);
		this->parameters->FindParam(&sealevTime,SmbSealevEnum,time);
	}
	else {
		TdiffTime=0;
		sealevTime=0;
	}

	/*Recover info at the vertices: */
	GetInputListOnVertices(&h[0],ThicknessEnum);
	GetInputListOnVertices(&s[0],SurfaceEnum);
	GetInputListOnVertices(&s0p[0],SmbS0pEnum);
	GetInputListOnVertices(&s0t[0],SmbS0tEnum);

	/*measure the surface mass balance*/
	for (int iv = 0; iv<numvertices; iv++){
		agd[iv]=PddSurfaceMassBalance(&monthlytemperatures[iv*12], &monthlyprec[iv*12],
					pdds, pds, &melt[iv], &accu[iv], signorm, yts, h[iv], s[iv],
					desfac, s0t[iv], s0p[iv],rlaps,rlapslgm,TdiffTime,sealevTime,
					rho_water,rho_ice);
		/*Get yearlytemperatures */
		for(int month=0;month<12;month++) {
			yearlytemperatures[iv]=yearlytemperatures[iv]+(monthlytemperatures[iv*12+month]+273.15)*mavg; // Has to be in Kelvin
		}
	}

	/*Update inputs*/
	// TransientInput* NewTemperatureInput = new TransientInput(SmbMonthlytemperaturesEnum);
	// TransientInput* NewPrecipitationInput = new TransientInput(SmbPrecipitationEnum);
	// for (int imonth=0;imonth<12;imonth++) {
	//   for(i=0;i<numvertices;i++) tmp[i]=monthlytemperatures[i*12+imonth];
	//   TriaInput* newmonthinput1 = new TriaInput(SmbMonthlytemperaturesEnum,&tmp[0],P1Enum);
	//   NewTemperatureInput->AddTimeInput(newmonthinput1,time+imonth/12.*yts);
	//
	//   for(i=0;i<numvertices;i++) tmp[i]=monthlyprec[i*12+imonth]/yts;
	//   TriaInput* newmonthinput2 = new TriaInput(SmbPrecipitationEnum,&tmp[0],P1Enum);
	//   NewPrecipitationInput->AddTimeInput(newmonthinput2,time+imonth/12.*yts);
	// }
	// NewTemperatureInput->Configure(this->parameters);
	// NewPrecipitationInput->Configure(this->parameters);

	switch(this->ObjectEnum()){
		case TriaEnum:  
			// this->inputs->AddInput(new TriaInput(TemperatureEnum,&yearlytemperatures[0],P1Enum));
			this->inputs->AddInput(new TriaInput(TemperaturePDDEnum,&yearlytemperatures[0],P1Enum));
			this->inputs->AddInput(new TriaInput(SmbMassBalanceEnum,&agd[0],P1Enum));
			this->inputs->AddInput(new TriaInput(SmbAccumulationEnum,&accu[0],P1Enum));
			this->inputs->AddInput(new TriaInput(SmbMeltEnum,&melt[0],P1Enum));
			break;
		case PentaEnum:
			if(IsOnSurface()){
				/*Here, we want to change the BC of the thermal model, keep
				 * the temperatures as they are for the base of the penta and
				 * yse yearlytemperatures for the top*/
				GetInputListOnVertices(&s[0],TemperatureEnum);
				yearlytemperatures[0] = s[0];
				yearlytemperatures[1] = s[1];
				yearlytemperatures[2] = s[2];
				this->inputs->AddInput(new PentaInput(TemperatureEnum,&yearlytemperatures[0],P1Enum));

				bool isenthalpy;
				this->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);
				if(isenthalpy){
					/*Convert that to enthalpy for the enthalpy model*/
					IssmDouble enthalpy[6];
					GetInputListOnVertices(&enthalpy[0],EnthalpyEnum);
					ThermalToEnthalpy(&enthalpy[3],yearlytemperatures[3],0.,0.);
					ThermalToEnthalpy(&enthalpy[4],yearlytemperatures[4],0.,0.);
					ThermalToEnthalpy(&enthalpy[5],yearlytemperatures[5],0.,0.);
					this->inputs->AddInput(new PentaInput(EnthalpyEnum,&enthalpy[0],P1Enum));
				}
			}
			this->inputs->AddInput(new PentaInput(SmbMassBalanceEnum,&agd[0],P1Enum));
			this->inputs->AddInput(new PentaInput(TemperaturePDDEnum,&yearlytemperatures[0],P1Enum));
			this->InputExtrude(TemperaturePDDEnum,-1);
			this->InputExtrude(SmbMassBalanceEnum,-1);
			break;
		case TetraEnum: 
			if(IsOnSurface()){
				GetInputListOnVertices(&s[0],TemperatureEnum);
				yearlytemperatures[0] = s[0];
				yearlytemperatures[1] = s[1];
				yearlytemperatures[2] = s[2];
				this->inputs->AddInput(new TetraInput(TemperatureEnum,&yearlytemperatures[0],P1Enum));
			}
			this->inputs->AddInput(new TetraInput(SmbMassBalanceEnum,&agd[0],P1Enum));
			this->inputs->AddInput(new TetraInput(TemperaturePDDEnum,&yearlytemperatures[0],P1Enum));
			this->InputExtrude(TemperaturePDDEnum,-1);
			this->InputExtrude(SmbMassBalanceEnum,-1);
			break;
		default: _error_("Not implemented yet");
	}
	// this->inputs->AddInput(NewTemperatureInput);
	// this->inputs->AddInput(NewPrecipitationInput);
	// this->inputs->AddInput(new TriaVertexInput(ThermalSpcTemperatureEnum,&Tsurf[0]));

	//this->InputExtrude(SmbMassBalanceEnum,-1);
	// this->InputExtrude(SmbMonthlytemperaturesEnum,-1);
	// this->InputExtrude(SmbPrecipitationEnum,-1);

	/*clean-up*/
	delete gauss;
	xDelete<IssmDouble>(monthlytemperatures);
	xDelete<IssmDouble>(monthlyprec);
	xDelete<IssmDouble>(agd);
	xDelete<IssmDouble>(melt);
	xDelete<IssmDouble>(accu);
	xDelete<IssmDouble>(yearlytemperatures);
	xDelete<IssmDouble>(h);
	xDelete<IssmDouble>(s);
	xDelete<IssmDouble>(s0t);
	xDelete<IssmDouble>(s0p);
	xDelete<IssmDouble>(tmp);

}
/*}}}*/
IssmDouble Element::PureIceEnthalpy(IssmDouble pressure){/*{{{*/
	return this->matpar->PureIceEnthalpy(pressure);
}/*}}}*/
void       Element::ResultInterpolation(int* pinterpolation,int* pnodesperelement,int* parray_size, int output_enum){/*{{{*/

	/*Some intputs need to be computed, even if they are already in inputs, they might not be up to date!*/
	switch(output_enum){
		case ViscousHeatingEnum: this->ViscousHeatingCreateInput(); break;
		case StressMaxPrincipalEnum: this->StressMaxPrincipalCreateInput(); break;
		case StressTensorxxEnum: 
		case StressTensorxyEnum: 
		case StressTensorxzEnum: 
		case StressTensoryyEnum: 
		case StressTensoryzEnum: 
		case StressTensorzzEnum: this->ComputeStressTensor(); break;
		case StrainRatexxEnum:
		case StrainRatexyEnum:
		case StrainRatexzEnum:
		case StrainRateyyEnum:
		case StrainRateyzEnum:
		case StrainRatezzEnum:
		case StrainRateeffectiveEnum: this->ComputeStrainRate(); break;
		case DeviatoricStressxxEnum: 
		case DeviatoricStressxyEnum: 
		case DeviatoricStressxzEnum: 
		case DeviatoricStressyyEnum: 
		case DeviatoricStressyzEnum: 
		case DeviatoricStresszzEnum: 
		case DeviatoricStresseffectiveEnum: this->ComputeDeviatoricStressTensor(); break;
		case EsaStrainratexxEnum:
		case EsaStrainratexyEnum: 
		case EsaStrainrateyyEnum: 
		case EsaRotationrateEnum: this->ComputeEsaStrainAndVorticity(); break;
		case SigmaNNEnum: this->ComputeSigmaNN(); break;
		case LambdaSEnum: this->ComputeLambdaS(); break;
		case NewDamageEnum: this->ComputeNewDamage(); break;
		case StressIntensityFactorEnum: this->StressIntensityFactor(); break;
		case CalvingratexEnum:
		case CalvingrateyEnum:
		case CalvingCalvingrateEnum:
			this->StrainRateparallel();
			this->StrainRateperpendicular();
			int calvinglaw;
			this->FindParam(&calvinglaw,CalvingLawEnum);
			switch(calvinglaw){
				case DefaultCalvingEnum:
					//do nothing
					break;
				case CalvingLevermannEnum:
					this->CalvingRateLevermann();
					break;
				case CalvingDevEnum:
					this->CalvingRateDev();
					break;
				default:
					_error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
			}
			break;
		case StrainRateparallelEnum: this->StrainRateparallel(); break;
		case StrainRateperpendicularEnum: this->StrainRateperpendicular(); break;
	}

	/*Find input*/
	Input* input=this->inputs->GetInput(output_enum);

	/*If this input is not already in Inputs, maybe it needs to be computed?*/
	if(!input) _error_("input "<<EnumToStringx(output_enum)<<" not found in element");

	/*Assign output pointer*/
	*pinterpolation   = input->GetResultInterpolation();
	*pnodesperelement = input->GetResultNumberOfNodes();
	*parray_size      = input->GetResultArraySize();
}/*}}}*/
void       Element::ResultToPatch(IssmDouble* values,int nodesperelement,int output_enum){/*{{{*/

	Input* input=this->inputs->GetInput(output_enum);
	if(!input) _error_("input "<<EnumToStringx(output_enum)<<" not found in element");

	input->ResultToPatch(values,nodesperelement,this->Sid());

} /*}}}*/
void       Element::ResultToMatrix(IssmDouble* values,int ncols,int output_enum){/*{{{*/

	Input* input=this->inputs->GetInput(output_enum);
	if(!input) _error_("input "<<EnumToStringx(output_enum)<<" not found in element");

	input->ResultToMatrix(values,ncols,this->Sid());

} /*}}}*/
void       Element::ResultToVector(Vector<IssmDouble>* vector,int output_enum){/*{{{*/

	Input* input=this->inputs->GetInput(output_enum);
	if(!input) _error_("input "<<EnumToStringx(output_enum)<<" not found in element");

	switch(input->GetResultInterpolation()){
		case P0Enum:{
			IssmDouble  value;
			bool        bvalue;
			Input*      input = this->GetInput(output_enum); _assert_(input);
			switch(input->ObjectEnum()){
				case DoubleInputEnum:
					input->GetInputValue(&value);
					break;
				case BoolInputEnum:
					input->GetInputValue(&bvalue);
					value=reCast<IssmDouble>(bvalue);
					break;
				default:
					Gauss* gauss = this->NewGauss();
					input->GetInputValue(&value,gauss);
					delete gauss;
			}
			vector->SetValue(this->Sid(),value,INS_VAL);
			break;
		}
		case P1Enum:{
			int         numvertices = this->GetNumberOfVertices();
			IssmDouble *values      = xNew<IssmDouble>(numvertices);
			int        *connectivity= xNew<int>(numvertices);
			int        *sidlist     = xNew<int>(numvertices);

			this->GetVerticesSidList(sidlist);
			this->GetVerticesConnectivityList(connectivity);
			this->GetInputListOnVertices(values,output_enum);
			for(int i=0;i<numvertices;i++) values[i] = values[i]/reCast<IssmDouble>(connectivity[i]);

			vector->SetValues(numvertices,sidlist,values,ADD_VAL);

			xDelete<IssmDouble>(values);
			xDelete<int>(connectivity);
			xDelete<int>(sidlist);
			break;
		}
		default:
					 _error_("interpolation "<<EnumToStringx(input->GetResultInterpolation())<<" not supported yet");
	}
} /*}}}*/
void       Element::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int set1_enum,int set2_enum){/*{{{*/

	/*Intermediaries*/
	const int numnodes = this->GetNumberOfNodes();

	/*Output */
	int d_nz = 0;
	int o_nz = 0;

	/*Loop over all nodes*/
	for(int i=0;i<numnodes;i++){

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

	/*Special case: 2d/3d coupling, the node of this element might be connected
	 *to the basal element*/
	int analysis_type,approximation,numlayers;
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	if(analysis_type==StressbalanceAnalysisEnum){
		inputs->GetInputValue(&approximation,ApproximationEnum);
		if(approximation==SSAHOApproximationEnum || approximation==SSAFSApproximationEnum){
			parameters->FindParam(&numlayers,MeshNumberoflayersEnum);
			o_nz += numlayers*3;
			d_nz += numlayers*3;
		}
	}

	/*Assign output pointers: */
	*pd_nz=d_nz;
	*po_nz=o_nz;
}
/*}}}*/
int        Element::Sid(){/*{{{*/

	return this->sid;

}
/*}}}*/
void       Element::SmbGemb(){/*{{{*/

	/*Intermediary variables: {{{*/
	IssmDouble isinitialized;
	IssmDouble zTop,dzTop,zMax,zMin,zY,dzMin;
	IssmDouble Tmean; 
	IssmDouble C; 
	IssmDouble Tz,Vz; 
	IssmDouble rho_ice, rho_water,aSnow,aIce;
	IssmDouble time,dt;
	IssmDouble t,smb_dt;
	IssmDouble yts;
	IssmDouble Ta,V,dlw,dsw,P,eAir,pAir;
	int        aIdx=0;
	int        denIdx=0;
	int        swIdx=0;
	IssmDouble cldFrac,t0wet, t0dry, K;
	IssmDouble ulw;
	IssmDouble netSW;
	IssmDouble netLW;
	IssmDouble lhf, shf, dayEC;
	IssmDouble initMass;
	IssmDouble sumR, sumM, sumEC, sumP, sumW,sumMassAdd;
	IssmDouble sumdz_add;
	IssmDouble sumMass,dMass;
	bool isgraingrowth,isalbedo,isshortwave,isthermal,isaccumulation,ismelt,isdensification,isturbulentflux;
	IssmDouble init_scaling;

	/*}}}*/
	/*Output variables:{{{ */
	IssmDouble* dz=NULL;
	IssmDouble* d = NULL;
	IssmDouble* re = NULL;
	IssmDouble* gdn = NULL;
	IssmDouble* gsp = NULL;
	IssmDouble  EC = 0;
	IssmDouble* W = NULL;
	IssmDouble* a = NULL;
	IssmDouble* swf=NULL;
	IssmDouble* T = NULL;
	IssmDouble  T_bottom;
	IssmDouble  M;
	IssmDouble  R; 
	IssmDouble  mAdd;
	IssmDouble  dz_add;
    
	IssmDouble* dzini=NULL;
	IssmDouble* dini = NULL;
	IssmDouble* reini = NULL;
	IssmDouble* gdnini = NULL;
	IssmDouble* gspini = NULL;
	IssmDouble* Wini = NULL;
	IssmDouble* aini = NULL;
	IssmDouble* Tini = NULL;
    
	int         m;
	int         count=0;
	/*}}}*/

	/*only compute SMB at the surface: */
	if (!IsOnSurface()) return;


	/*Retrieve material properties and parameters:{{{ */
	rho_ice = matpar->GetMaterialParameter(MaterialsRhoIceEnum);
	rho_water = matpar->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	parameters->FindParam(&aSnow,SmbASnowEnum);
	parameters->FindParam(&aIce,SmbAIceEnum);
	parameters->FindParam(&time,TimeEnum);                        /*transient core time at which we run the smb core*/
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);          /*transient core time step*/
	parameters->FindParam(&yts,ConstantsYtsEnum);
	parameters->FindParam(&smb_dt,SmbDtEnum);                     /*time period for the smb solution,  usually smaller than the glaciological dt*/
	parameters->FindParam(&aIdx,SmbAIdxEnum);
	parameters->FindParam(&denIdx,SmbDenIdxEnum);
	parameters->FindParam(&swIdx,SmbSwIdxEnum);
	parameters->FindParam(&cldFrac,SmbCldFracEnum);
	parameters->FindParam(&t0wet,SmbT0wetEnum);
	parameters->FindParam(&t0dry,SmbT0dryEnum);
	parameters->FindParam(&K,SmbKEnum);
	parameters->FindParam(&isgraingrowth,SmbIsgraingrowthEnum);
	parameters->FindParam(&isalbedo,SmbIsalbedoEnum);
	parameters->FindParam(&isshortwave,SmbIsshortwaveEnum);
	parameters->FindParam(&isthermal,SmbIsthermalEnum);
	parameters->FindParam(&isaccumulation,SmbIsaccumulationEnum);
	parameters->FindParam(&ismelt,SmbIsmeltEnum);
	parameters->FindParam(&isdensification,SmbIsdensificationEnum);
	parameters->FindParam(&isturbulentflux,SmbIsturbulentfluxEnum);
	parameters->FindParam(&init_scaling,SmbInitDensityScalingEnum);
    
	/*}}}*/
	/*Retrieve inputs: {{{*/
	Input* zTop_input=this->GetInput(SmbZTopEnum); _assert_(zTop_input); 
	Input* dzTop_input=this->GetInput(SmbDzTopEnum); _assert_(dzTop_input); 
	Input* dzMin_input=this->GetInput(SmbDzMinEnum); _assert_(dzMin_input); 
	Input* zMax_input=this->GetInput(SmbZMaxEnum); _assert_(zMax_input); 
	Input* zMin_input=this->GetInput(SmbZMinEnum); _assert_(zMin_input); 
	Input* zY_input=this->GetInput(SmbZYEnum); _assert_(zY_input); 
	Input* Tmean_input=this->GetInput(SmbTmeanEnum); _assert_(Tmean_input);
	Input* C_input=this->GetInput(SmbCEnum); _assert_(C_input);
	Input* Tz_input=this->GetInput(SmbTzEnum); _assert_(Tz_input);
	Input* Vz_input=this->GetInput(SmbVzEnum); _assert_(Vz_input);
	Input* Ta_input=this->GetInput(SmbTaEnum); _assert_(Ta_input);
	Input* V_input=this->GetInput(SmbVEnum); _assert_(V_input);
	Input* Dlwr_input=this->GetInput(SmbDlwrfEnum); _assert_(Dlwr_input);
	Input* Dswr_input=this->GetInput(SmbDswrfEnum); _assert_(Dswr_input);
	Input* P_input=this->GetInput(SmbPEnum); _assert_(P_input);
	Input* eAir_input=this->GetInput(SmbEAirEnum); _assert_(eAir_input);
	Input* pAir_input=this->GetInput(SmbPAirEnum); _assert_(pAir_input);
	Input* isinitialized_input=this->GetInput(SmbIsInitializedEnum); _assert_(isinitialized_input);
	/*Retrieve input values:*/
	Gauss* gauss=this->NewGauss(1); gauss->GaussPoint(0);

	zTop_input->GetInputValue(&zTop,gauss);
	dzTop_input->GetInputValue(&dzTop,gauss);
	dzMin_input->GetInputValue(&dzMin,gauss);
	zMax_input->GetInputValue(&zMax,gauss); 
	zMin_input->GetInputValue(&zMin,gauss); 
	zY_input->GetInputValue(&zY,gauss);
	Tmean_input->GetInputValue(&Tmean,gauss);
	C_input->GetInputValue(&C,gauss);
	Tz_input->GetInputValue(&Tz,gauss);
	Vz_input->GetInputValue(&Vz,gauss);
	isinitialized_input->GetInputValue(&isinitialized);
	/*}}}*/

	/*First, check that the initial structures have been setup in GEMB. If not, initialize profile variables: layer thickness dz, * density d, temperature T, etc. {{{*/
	if(isinitialized==0.0){
        if(VerboseSmb() && this->Sid()==0)_printf0_("smb core: Initializing grid\n");
        //if(this->Sid()==1) for(int i=0;i<m;i++)_printf_("z[" << i << "]=" <<
        //dz[i] << "\n");
        
        DoubleArrayInput* dz_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbDziniEnum)); _assert_(dz_input);
        DoubleArrayInput* d_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbDiniEnum));_assert_(d_input);
        DoubleArrayInput* re_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbReiniEnum));_assert_(re_input);
        DoubleArrayInput* gdn_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbGdniniEnum));_assert_(gdn_input);
        DoubleArrayInput* gsp_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbGspiniEnum));_assert_(gsp_input);
        DoubleInput* EC_input= dynamic_cast<DoubleInput*>(this->GetInput(SmbECiniEnum));_assert_(EC_input);
        DoubleArrayInput* W_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbWiniEnum));_assert_(W_input);
        DoubleArrayInput* a_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbAiniEnum));_assert_(a_input);
        DoubleArrayInput* T_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbTiniEnum));_assert_(T_input);

        dz_input->GetValues(&dzini,&m);
        d_input->GetValues(&dini,&m);
        re_input->GetValues(&reini,&m);
        gdn_input->GetValues(&gdnini,&m);
        gsp_input->GetValues(&gspini,&m);
        EC_input->GetInputValue(&EC);
        W_input->GetValues(&Wini,&m);
        a_input->GetValues(&aini,&m);
        T_input->GetValues(&Tini,&m);
        
        /*Retrive the correct value of m (without the zeroes at the end)*/
        Input* Size_input=this->GetInput(SmbSizeiniEnum); _assert_(Size_input);
        Size_input->GetInputValue(&m);
        
        if(m==2){ //Snow properties are initialized with default values. Vertical grid has to be initialized too
//            if(VerboseSmb() && this->Sid()==0)_printf0_("Snow properties initialized w DEFAULT values\n");
            
            /*initialize profile variables:*/
            GembgridInitialize(&dz, &m, zTop, dzTop, zMax, zY);
            
            d = xNewZeroInit<IssmDouble>(m); for(int i=0;i<m;i++)d[i]=dini[0]; //ice density [kg m-3]
            re = xNewZeroInit<IssmDouble>(m); for(int i=0;i<m;i++)re[i]=reini[0];         //set grain size to old snow [mm]
            gdn = xNewZeroInit<IssmDouble>(m); for(int i=0;i<m;i++)gdn[i]=gdnini[0];         //set grain dentricity to old snow
            gsp = xNewZeroInit<IssmDouble>(m); for(int i=0;i<m;i++)gsp[i]=gspini[0];         //set grain sphericity to old snow
            W = xNewZeroInit<IssmDouble>(m); for(int i=0;i<m;i++)W[i]=Wini[0];             //set water content to zero [kg m-2]
            a = xNewZeroInit<IssmDouble>(m); for(int i=0;i<m;i++)a[i]=aini[0];         //set albedo equal to fresh snow [fraction]
            T = xNewZeroInit<IssmDouble>(m); for(int i=0;i<m;i++)T[i]=Tmean;         //set initial grid cell temperature to the annual mean temperature [K]
/*            /!\ Default value of T can not be retrived from SMBgemb.m (like other snow properties) because don't know Tmean yet when set default values.
            Default value of 0C given in SMBgemb.m is overwritten here with value of Tmean*/
            
            //fixed lower temperature bounday condition - T is fixed
            T_bottom=T[m-1];
        }
        else{ //Retrieve snow properties from previous run. Need to provide values for all layers
//            if(VerboseSmb() && this->Sid()==0)_printf0_("Snow properties initialized w RESTART values\n");
            
            dz = xNewZeroInit<IssmDouble>(m);for(int i=0;i<m;i++)dz[i]=dzini[i];
            d = xNewZeroInit<IssmDouble>(m);for(int i=0;i<m;i++)d[i]=dini[i];
            re = xNewZeroInit<IssmDouble>(m);for(int i=0;i<m;i++)re[i]=reini[i];
            gdn = xNewZeroInit<IssmDouble>(m);for(int i=0;i<m;i++)gdn[i]=gdnini[i];
            gsp = xNewZeroInit<IssmDouble>(m);for(int i=0;i<m;i++)gsp[i]=gspini[i];
            W = xNewZeroInit<IssmDouble>(m);for(int i=0;i<m;i++)W[i]=Wini[i];
            a = xNewZeroInit<IssmDouble>(m);for(int i=0;i<m;i++)a[i]=aini[i];
            T = xNewZeroInit<IssmDouble>(m);for(int i=0;i<m;i++)T[i]=Tini[i];

            //fixed lower temperature bounday condition - T is fixed
            T_bottom=T[m-1];
        }
        
        /*Flag the initialization:*/
        this->AddInput(new DoubleInput(SmbIsInitializedEnum,1.0));
    }
    else{
        /*Recover inputs: */
        DoubleArrayInput* dz_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbDzEnum)); _assert_(dz_input);
        DoubleArrayInput* d_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbDEnum));_assert_(d_input);
        DoubleArrayInput* re_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbReEnum));_assert_(re_input);
        DoubleArrayInput* gdn_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbGdnEnum));_assert_(gdn_input);
        DoubleArrayInput* gsp_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbGspEnum));_assert_(gsp_input);
        DoubleInput* EC_input= dynamic_cast<DoubleInput*>(this->GetInput(SmbECEnum));_assert_(EC_input);
        DoubleArrayInput* W_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbWEnum));_assert_(W_input);
        DoubleArrayInput* a_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbAEnum));_assert_(a_input);
        DoubleArrayInput* T_input= dynamic_cast<DoubleArrayInput*>(this->GetInput(SmbTEnum));_assert_(T_input);

        /*Recover arrays: */
        dz_input->GetValues(&dz,&m);
        d_input->GetValues(&d,&m);
        re_input->GetValues(&re,&m);
        gdn_input->GetValues(&gdn,&m);
        gsp_input->GetValues(&gsp,&m);
        EC_input->GetInputValue(&EC);
        W_input->GetValues(&W,&m);
        a_input->GetValues(&a,&m);
        T_input->GetValues(&T,&m);

        //fixed lower temperature bounday condition - T is fixed
        T_bottom=T[m-1];

    } /*}}}*/

	// determine initial mass [kg]
	initMass=0; for(int i=0;i<m;i++) initMass += dz[i]*d[i] + W[i];
    
        // initialize cumulative variables
	sumR = 0; sumM = 0; sumEC = 0; sumP = 0; sumMassAdd = 0;
	sumdz_add=0;

	//before starting loop, realize that the transient core runs this smb_core at time = time +deltaT. 
	//go back to time - deltaT: 
	time-=dt;

	/*Start loop: */
	count=1;
	for (t=time;t<time+dt;t=t+smb_dt){

		if(VerboseSmb() && this->Sid()==0 && IssmComm::GetRank()==0)_printf0_("Time: t=" << setprecision(8) << t/365.0/24.0/3600.0 << " yr/" << (time+dt)/365.0/24.0/3600.0 << " yr" << setprecision(3) << " Step: " << count << "\n");

		/*extract daily data:{{{*/
		Ta_input->GetInputValue(&Ta,gauss,t);//screen level air temperature [K]
		V_input->GetInputValue(&V,gauss,t);  //wind speed [m s-1]
		Dlwr_input->GetInputValue(&dlw,gauss,t);   //downward longwave radiation flux [W m-2]
		Dswr_input->GetInputValue(&dsw,gauss,t);   //downward shortwave radiation flux [W m-2]
		P_input->GetInputValue(&P,gauss,t);        //precipitation [kg m-2]
		eAir_input->GetInputValue(&eAir,gauss,t);  //screen level vapor pressure [Pa]
		pAir_input->GetInputValue(&pAir,gauss,t);  // screen level air pressure [Pa]
		//_printf_("Time: " << t << " Ta: " << Ta << " V: " << V << " dlw: " << dlw << " dsw: " << dsw << " P: " << P << " eAir: " << eAir << " pAir: " << pAir << "\n");
		/*}}}*/

		/*Snow grain metamorphism:*/
		if(isgraingrowth)grainGrowth(re, gdn, gsp, T, dz, d, W, smb_dt, m, aIdx,this->Sid());

		/*Snow, firn and ice albedo:*/
		if(isalbedo)albedo(a,aIdx,re,d,cldFrac,aIce, aSnow,T,W,P,EC,t0wet,t0dry,K,smb_dt,m,this->Sid());
		
					
		/*Distribution of absorbed short wave radation with depth:*/
		if(isshortwave)shortwave(&swf, swIdx, aIdx, dsw, a[0], d, dz, re,m,this->Sid());
		
		/*Calculate net shortwave [W m-2]*/
		netSW = cellsum(swf,m);

		/*Thermal profile computation:*/
		if(isthermal)thermo(&EC, T, dz, d, swf, dlw, Ta, V, eAir, pAir, W[0], smb_dt, m, Vz, Tz,this->Sid());

		/*Change in thickness of top cell due to evaporation/condensation  assuming same density as top cell. 
		 * need to fix this in case all or more of cell evaporates */
		dz[0] = dz[0] + EC / d[0];
		
		/*Add snow/rain to top grid cell adjusting cell depth, temperature and density*/
		if(isaccumulation)accumulation(&T, &dz, &d, &W, &a, &re, &gdn, &gsp, &m, Ta, P, dzMin, aSnow,this->Sid());

		/*Calculate water production, M [kg m-2] resulting from snow/ice temperature exceeding 273.15 deg K
		 * (> 0 deg C), runoff R [kg m-2] and resulting changes in density and determine wet compaction [m]*/
		if(ismelt)melt(&M, &R, &mAdd, &dz_add, &T, &d, &dz, &W, &a, &re, &gdn, &gsp, &m, dzMin, zMax, zMin, zTop,this->Sid());

		/*Allow non-melt densification and determine compaction [m]*/
		if(isdensification)densification(d,dz, T, re, denIdx, C, smb_dt, Tmean,rho_ice,m,this->Sid());
		
		/*Calculate upward longwave radiation flux [W m-2] not used in energy balance. Calculated for every 
		 * sub-time step in thermo equations*/
		ulw = 5.67E-8 * pow(T[0],4.0);

		/*Calculate net longwave [W m-2]*/
		netLW = dlw - ulw;
		
		/*Calculate turbulent heat fluxes [W m-2]*/
		if(isturbulentflux)turbulentFlux(&shf, &lhf, &dayEC, Ta, T[0], V, eAir, pAir, d[0], W[0], Vz, Tz,this->Sid());
		
		/*Verbose some resuls in debug mode: {{{*/
		if(VerboseSmb() && 0){ 
			_printf_("smb log: count[" << count << "] m[" << m << "] " 
				<< setprecision(16)   << "T[" << cellsum(T,m)  << "] " 
					                  << "d[" << cellsum(d,m)  << "] "
					                  << "dz[" << cellsum(dz,m)  << "] "
					                  << "a[" << cellsum(a,m)  << "] "
					                  << "W[" << cellsum(W,m)  << "] "
					                  << "re[" << cellsum(re,m)  << "] "
					                  << "gdn[" << cellsum(gdn,m)  << "] "
					                  << "gsp[" << cellsum(gsp,m)  << "] "
					                  << "swf[" << netSW << "] "
									  << "\n");
		} /*}}}*/
		
		/*Sum component mass changes [kg m-2]*/
		sumMassAdd = mAdd + sumMassAdd;
		sumM = M + sumM;
		sumR = R + sumR;
		sumW = cellsum(W,m);
		sumP = P +  sumP;
		sumEC = sumEC + EC;  // evap (-)/cond(+)
		sumdz_add=dz_add+sumdz_add;

		/*Calculate total system mass:*/
		sumMass=0; for(int i=0;i<m;i++) sumMass += dz[i]*d[i];

		#ifndef _HAVE_ADOLC_ //we want to avoid the round operation at all cost. Not differentiable.
		dMass = sumMass + sumR + sumW - sumP - sumEC - initMass - sumMassAdd;
		dMass = round(dMass * 100.0)/100.0;

		/*Check mass conservation:*/
		if (dMass != 0.0) _printf_("total system mass not conserved in MB function \n");
		#endif
		
		/*Check bottom grid cell T is unchanged:*/
		if (T[m-1]!=T_bottom) _printf_("T(end)~=T_bottom" << "\n");
		
		/*Free ressources: */
		xDelete<IssmDouble>(swf);

		/*increase counter:*/
		count++;
	} //for (t=time;t<time+dt;t=t+smb_dt)

	/*Save generated inputs: */
	this->AddInput(new DoubleArrayInput(SmbDzEnum,dz,m));
	this->AddInput(new DoubleArrayInput(SmbDEnum,d,m));
	this->AddInput(new DoubleArrayInput(SmbReEnum,re,m));
	this->AddInput(new DoubleArrayInput(SmbGdnEnum,gdn,m));
	this->AddInput(new DoubleArrayInput(SmbGspEnum,gsp,m));
	this->AddInput(new DoubleArrayInput(SmbTEnum,T,m));
	this->AddInput(new DoubleInput(SmbECEnum,sumEC/yts));
	this->AddInput(new DoubleArrayInput(SmbWEnum,W,m));
	this->AddInput(new DoubleArrayInput(SmbAEnum,a,m));
	this->AddInput(new DoubleInput(SmbMassBalanceEnum,(sumP + sumEC -sumR)/yts));
	this->AddInput(new DoubleInput(SmbRunoffEnum,sumR/yts));
	this->AddInput(new DoubleInput(SmbPrecipitationEnum,sumP/yts));
	this->AddInput(new DoubleInput(SmbDz_addEnum,sumdz_add/yts));
	this->AddInput(new DoubleInput(SmbM_addEnum,sumMassAdd/yts));

	/*Free allocations:{{{*/
	xDelete<IssmDouble>(dz);
	xDelete<IssmDouble>(d);
	xDelete<IssmDouble>(re);
	xDelete<IssmDouble>(gdn);
	xDelete<IssmDouble>(gsp);
	xDelete<IssmDouble>(W);
	xDelete<IssmDouble>(a);
	xDelete<IssmDouble>(T);
	xDelete<IssmDouble>(dzini);
	xDelete<IssmDouble>(dini);
	xDelete<IssmDouble>(reini);
	xDelete<IssmDouble>(gdnini);
	xDelete<IssmDouble>(gspini);
	xDelete<IssmDouble>(Wini);
	xDelete<IssmDouble>(aini);
	xDelete<IssmDouble>(Tini);

	delete gauss;
	/*}}}*/
}
/*}}}*/
void       Element::StrainRateESA(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble dvx[3];
	IssmDouble dvy[3];

	/*Check that both inputs have been found*/
	if(!vx_input || !vy_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input << ", vy: " << vy_input << "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	epsilon[0] = dvx[0];	// normal strain rate x-direction  
	epsilon[1] = dvy[1]; // normal strain rate y-direction 
	epsilon[2] = 0.5*(dvx[1] + dvy[0]); // shear strain rate 
	epsilon[3] = 0.5*(dvx[1] - dvy[0]); // rotation rate 

}/*}}}*/
void       Element::StrainRateFS(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/
	/*Compute the 3d Strain Rate (6 components):
	 *
	 * epsilon=[exx eyy ezz exy exz eyz]
	 */

	/*Intermediaries*/
	IssmDouble dvx[3];
	IssmDouble dvy[3];
	IssmDouble dvz[3];

	/*Check that both inputs have been found*/
	if (!vx_input || !vy_input || !vz_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input << ", vy: " << vy_input << ", vz: " << vz_input << "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
	epsilon[0] = dvx[0];
	epsilon[1] = dvy[1];
	epsilon[2] = dvz[2];
	epsilon[3] = 0.5*(dvx[1] + dvy[0]);
	epsilon[4] = 0.5*(dvx[2] + dvz[0]);
	epsilon[5] = 0.5*(dvy[2] + dvz[1]);

}/*}}}*/
void       Element::StrainRateHO(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/
	/*Compute the 3d Blatter/HOStrain Rate (5 components):
	 *
	 * epsilon=[exx eyy exy exz eyz]
	 *
	 * with exz=1/2 du/dz
	 *      eyz=1/2 dv/dz
	 *
	 * the contribution of vz is neglected
	 */

	/*Intermediaries*/
	IssmDouble dvx[3];
	IssmDouble dvy[3];

	/*Check that both inputs have been found*/
	if (!vx_input || !vy_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input << ", vy: " << vy_input << "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	epsilon[0] = dvx[0];
	epsilon[1] = dvy[1];
	epsilon[2] = 0.5*(dvx[1] + dvy[0]);
	epsilon[3] = 0.5*dvx[2];
	epsilon[4] = 0.5*dvy[2];

}/*}}}*/
void       Element::StrainRateHO2dvertical(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/
	/*Compute the 2d Blatter/HOStrain Rate (2 components):
	 *
	 * epsilon=[exx exz]
	 *
	 * with exz=1/2 du/dz
	 *
	 * the contribution of vz is neglected
	 */

	/*Intermediaries*/
	IssmDouble dvx[3];

	/*Check that both inputs have been found*/
	if (!vx_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input <<"\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	epsilon[0] = dvx[0];
	epsilon[1] = 0.5*dvx[1];

}/*}}}*/
void       Element::StrainRateSSA(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble dvx[3];
	IssmDouble dvy[3];

	/*Check that both inputs have been found*/
	if(!vx_input || !vy_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input << ", vy: " << vy_input << "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
	epsilon[0] = dvx[0];
	epsilon[1] = dvy[1];
	epsilon[2] = 0.5*(dvx[1] + dvy[0]);

}/*}}}*/
void       Element::StrainRateSSA1d(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble dvx[3];

	/*Check that both inputs have been found*/
	if (!vx_input){
		_error_("Input missing. Here are the input pointers we have for vx: " << vx_input << "\n");
	}

	/*Get strain rate assuming that epsilon has been allocated*/
	vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
	*epsilon = dvx[0];

}/*}}}*/
void       Element::StressMaxPrincipalCreateInput(void){/*{{{*/

	/*Intermediaries*/
	IssmDouble *xyz_list = NULL;
	IssmDouble  sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
	IssmDouble  a,b,c,d,x[3],max;
	int         dim,numroots;

	/*First: get stress tensor*/
	this->ComputeStressTensor();

	/*Get domain dimension*/
	this->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number vertices and allocate memory*/
	int         numvertices  = this->GetNumberOfVertices();
	IssmDouble* maxprincipal = xNew<IssmDouble>(numvertices);

	/*Retrieve all inputs and parameters*/
	this->GetVerticesCoordinatesBase(&xyz_list);
	Input* sigma_xx_input  = this->GetInput(StressTensorxxEnum); _assert_(sigma_xx_input);
	Input* sigma_yy_input  = this->GetInput(StressTensoryyEnum); _assert_(sigma_yy_input);
	Input* sigma_xy_input  = this->GetInput(StressTensorxyEnum); _assert_(sigma_xy_input);
	Input* sigma_xz_input  = NULL;
	Input* sigma_yz_input  = NULL;
	Input* sigma_zz_input  = NULL;
	if(dim==3){
		sigma_xz_input  = this->GetInput(StressTensorxzEnum); _assert_(sigma_xz_input);
		sigma_yz_input  = this->GetInput(StressTensoryzEnum); _assert_(sigma_yz_input);
		sigma_zz_input  = this->GetInput(StressTensorzzEnum); _assert_(sigma_zz_input);
	}

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for (int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

		sigma_xx_input->GetInputValue(&sigma_xx,gauss);
		sigma_yy_input->GetInputValue(&sigma_yy,gauss);
		sigma_xy_input->GetInputValue(&sigma_xy,gauss);
		if(dim==3){
			sigma_xz_input->GetInputValue(&sigma_xz,gauss);
			sigma_yz_input->GetInputValue(&sigma_yz,gauss);
			sigma_zz_input->GetInputValue(&sigma_zz,gauss);
		}

		if(dim==2){
			a = 0.;
			b = 1.;
			c = -sigma_yy -sigma_xx;
			d = sigma_xx*sigma_yy - sigma_xy*sigma_xy;
		}
		else{
			a = -1.;
			b = sigma_xx+sigma_yy+sigma_zz;
			c = -sigma_xx*sigma_yy -sigma_xx*sigma_zz -sigma_yy*sigma_zz + sigma_xy*sigma_xy +sigma_xz*sigma_xz +sigma_yz*sigma_yz;
			d = sigma_xx*sigma_yy*sigma_zz - sigma_xx*sigma_yz*sigma_yz -sigma_yy*sigma_xz*sigma_xz - sigma_zz*sigma_xy*sigma_xy + 2.*sigma_xy*sigma_xz*sigma_yz;
		}

		/*Get roots of polynomials*/
		cubic(a,b,c,d,x,&numroots);

		/*Initialize maximum eigne value*/
		if(numroots>0){
			max = fabs(x[0]);
		}
		else{
			_error_("No eigen value found");
		}

		/*Get max*/
		for(int i=1;i<numroots;i++){
			if(fabs(x[i])>max) max = fabs(x[i]);
		}

		maxprincipal[iv]=max;
	}

	/*Create input*/
	this->AddInput(StressMaxPrincipalEnum,maxprincipal,P1Enum);

	/*Clean up and return*/
	xDelete<IssmDouble>(maxprincipal);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
}
/*}}}*/
void       Element::ThermalToEnthalpy(IssmDouble* penthalpy,IssmDouble temperature,IssmDouble waterfraction,IssmDouble pressure){/*{{{*/
	matpar->ThermalToEnthalpy(penthalpy,temperature,waterfraction,pressure);
}/*}}}*/
IssmDouble Element::TMeltingPoint(IssmDouble pressure){/*{{{*/
	_assert_(matpar);
	return this->matpar->TMeltingPoint(pressure);
}/*}}}*/
void       Element::TransformInvStiffnessMatrixCoord(ElementMatrix* Ke,int transformenum){/*{{{*/

	/*All nodes have the same Coordinate System*/
	int numnodes  = this->GetNumberOfNodes();
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;

	/*Call core*/
	TransformInvStiffnessMatrixCoord(Ke,this->nodes,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformInvStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes_list,int numnodes,int* cs_array){/*{{{*/

	int         i,j;
	int         numdofs   = 0;
	IssmDouble *transform = NULL;
	IssmDouble *values    = NULL;

	/*Get total number of dofs*/
	for(i=0;i<numnodes;i++){
		switch(cs_array[i]){
			case PressureEnum: numdofs+=1; break;
			case XYEnum:       numdofs+=2; break;
			case XYZEnum:      numdofs+=3; break;
			default: _error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Copy current stiffness matrix*/
	values=xNew<IssmDouble>(Ke->nrows*Ke->ncols);
	for(i=0;i<Ke->nrows;i++) for(j=0;j<Ke->ncols;j++) values[i*Ke->ncols+j]=Ke->values[i*Ke->ncols+j];

	/*Get Coordinate Systems transform matrix*/
	CoordinateSystemTransform(&transform,nodes_list,numnodes,cs_array);

	/*Transform matrix: R*Ke*R^T */
	TripleMultiply(transform,numdofs,numdofs,0,
				values,Ke->nrows,Ke->ncols,0,
				transform,numdofs,numdofs,1,
				&Ke->values[0],0);

	/*Free Matrix*/
	xDelete<IssmDouble>(transform);
	xDelete<IssmDouble>(values);
}/*}}}*/
void       Element::TransformLoadVectorCoord(ElementVector* pe,int transformenum){/*{{{*/

	/*All nodes have the same Coordinate System*/
	int  numnodes = this->GetNumberOfNodes();
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;

	/*Call core*/
	this->TransformLoadVectorCoord(pe,this->nodes,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformLoadVectorCoord(ElementVector* pe,int* cs_array){/*{{{*/

	this->TransformLoadVectorCoord(pe,this->nodes,this->GetNumberOfNodes(),cs_array);

}/*}}}*/
void       Element::TransformLoadVectorCoord(ElementVector* pe,Node** nodes_list,int numnodes,int* cs_array){/*{{{*/

	int         i;
	int         numdofs   = 0;
	IssmDouble *transform = NULL;
	IssmDouble *values    = NULL;

	/*Get total number of dofs*/
	for(i=0;i<numnodes;i++){
		switch(cs_array[i]){
			case PressureEnum: numdofs+=1; break;
			case XYEnum:       numdofs+=2; break;
			case XYZEnum:      numdofs+=3; break;
			default: _error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Copy current load vector*/
	values=xNew<IssmDouble>(pe->nrows);
	for(i=0;i<pe->nrows;i++) values[i]=pe->values[i];

	/*Get Coordinate Systems transform matrix*/
	CoordinateSystemTransform(&transform,nodes_list,numnodes,cs_array);

	/*Transform matrix: R^T*pe */
	MatrixMultiply(transform,numdofs,numdofs,1,
				values,pe->nrows,1,0,
				&pe->values[0],0);

	/*Free Matrices*/
	xDelete<IssmDouble>(transform);
	xDelete<IssmDouble>(values);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* values,int transformenum){/*{{{*/

	/*All nodes have the same Coordinate System*/
	int  numnodes = this->GetNumberOfNodes();
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;

	/*Call core*/
	this->TransformSolutionCoord(values,this->nodes,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* values,int* transformenum_list){/*{{{*/
	this->TransformSolutionCoord(values,this->nodes,this->GetNumberOfNodes(),transformenum_list);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* values,int numnodes,int transformenum){/*{{{*/

	/*All nodes have the same Coordinate System*/
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;

	/*Call core*/
	this->TransformSolutionCoord(values,this->nodes,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* solution,int numnodes,int* cs_array){/*{{{*/
	this->TransformSolutionCoord(solution,this->nodes,numnodes,cs_array);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* values,Node** nodes_list,int numnodes,int transformenum){/*{{{*/
	/*NOT NEEDED*/
	/*All nodes have the same Coordinate System*/
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;

	/*Call core*/
	this->TransformSolutionCoord(values,nodes_list,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformSolutionCoord(IssmDouble* solution,Node** nodes_list,int numnodes,int* cs_array){/*{{{*/

	int         i;
	int         numdofs   = 0;
	IssmDouble *transform = NULL;
	IssmDouble *values    = NULL;

	/*Get total number of dofs*/
	for(i=0;i<numnodes;i++){
		switch(cs_array[i]){
			case PressureEnum: numdofs+=1; break;
			case XYEnum:       numdofs+=2; break;
			case XYZEnum:      numdofs+=3; break;
			default: _error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Copy current solution vector*/
	values=xNew<IssmDouble>(numdofs);
	for(i=0;i<numdofs;i++) values[i]=solution[i];

	/*Get Coordinate Systems transform matrix*/
	CoordinateSystemTransform(&transform,nodes_list,numnodes,cs_array);

	/*Transform matrix: R*U */
	MatrixMultiply(transform,numdofs,numdofs,0,
				values,numdofs,1,0,
				&solution[0],0);

	/*Free Matrices*/
	xDelete<IssmDouble>(transform);
	xDelete<IssmDouble>(values);
}/*}}}*/
void       Element::TransformStiffnessMatrixCoord(ElementMatrix* Ke,int transformenum){/*{{{*/

	/*All nodes have the same Coordinate System*/
	int  numnodes = this->GetNumberOfNodes();
	int* cs_array = xNew<int>(numnodes);
	for(int i=0;i<numnodes;i++) cs_array[i]=transformenum;

	/*Call core*/
	this->TransformStiffnessMatrixCoord(Ke,this->nodes,numnodes,cs_array);

	/*Clean-up*/
	xDelete<int>(cs_array);
}/*}}}*/
void       Element::TransformStiffnessMatrixCoord(ElementMatrix* Ke,int* transformenum_list){/*{{{*/
	this->TransformStiffnessMatrixCoord(Ke,this->nodes,this->GetNumberOfNodes(),transformenum_list);
}/*}}}*/
void       Element::TransformStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes_list,int numnodes,int* cs_array){/*{{{*/

	int         numdofs = 0;
	IssmDouble *transform = NULL;
	IssmDouble *values    = NULL;

	/*Get total number of dofs*/
	for(int i=0;i<numnodes;i++){
		switch(cs_array[i]){
			case PressureEnum: numdofs+=1; break;
			case XYEnum:       numdofs+=2; break;
			case XYZEnum:      numdofs+=3; break;
			default: _error_("Coordinate system " << EnumToStringx(cs_array[i]) << " not supported yet");
		}
	}

	/*Copy current stiffness matrix*/
	values=xNew<IssmDouble>(Ke->nrows*Ke->ncols);
	for(int i=0;i<Ke->nrows*Ke->ncols;i++) values[i]=Ke->values[i];

	/*Get Coordinate Systems transform matrix*/
	CoordinateSystemTransform(&transform,nodes_list,numnodes,cs_array);

	/*Transform matrix: R^T*Ke*R */
	TripleMultiply(transform,numdofs,numdofs,1,
				values,Ke->nrows,Ke->ncols,0,
				transform,numdofs,numdofs,0,
				&Ke->values[0],0);

	/*Free Matrix*/
	xDelete<IssmDouble>(transform);
	xDelete<IssmDouble>(values);
}/*}}}*/
void       Element::ViscousHeatingCreateInput(void){/*{{{*/

	/*Intermediaries*/
	IssmDouble phi;
	IssmDouble viscosity;
	IssmDouble epsilon[6];
	IssmDouble thickness;
	IssmDouble *xyz_list = NULL;

	/*Fetch number vertices and allocate memory*/
	int         numvertices    = this->GetNumberOfVertices();
	IssmDouble* viscousheating = xNew<IssmDouble>(numvertices);

	/*Retrieve all inputs and parameters*/
	this->GetVerticesCoordinatesBase(&xyz_list);
	Input* vx_input        = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input        = this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input        = this->GetInput(VzEnum); _assert_(vz_input);
	Input* thickness_input = this->GetInput(ThicknessEnum); _assert_(thickness_input);

	/*loop over vertices: */
	Gauss* gauss=this->NewGauss();
	for (int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

		thickness_input->GetInputValue(&thickness,gauss);

		this->StrainRateFS(&epsilon[0],xyz_list,gauss,vx_input,vy_input,vz_input);
		this->material->ViscosityFS(&viscosity,3,xyz_list,gauss,vx_input,vy_input,vz_input);
		this->GetPhi(&phi,&epsilon[0],viscosity);

		viscousheating[iv]=phi*thickness;
	}

	/*Create PentaVertex input, which will hold the basal friction:*/
	this->AddInput(ViscousHeatingEnum,viscousheating,P1Enum);

	/*Clean up and return*/
	xDelete<IssmDouble>(viscousheating);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
}
/*}}}*/
