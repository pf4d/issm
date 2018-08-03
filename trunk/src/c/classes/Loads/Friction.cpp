/*!\file Friction.c
 * \brief: implementation of the Friction object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "shared/shared.h"
/*}}}*/	

/*Constructors/destructors*/
Friction::Friction(){/*{{{*/
	this->element=NULL;
	this->dim=0;
	this->law=0;

}
/*}}}*/
Friction::Friction(Element* element_in,int dim_in){/*{{{*/

	this->element=element_in;
	this->dim=dim_in;
	element_in->FindParam(&this->law,FrictionLawEnum);
}
/*}}}*/
Friction::~Friction(){/*{{{*/
}
/*}}}*/

/*methods: */
void Friction::Echo(void){/*{{{*/
	_printf_("Friction:\n");
	_printf_("   dim: " << this->dim<< "\n");
}
/*}}}*/
void Friction::GetAlphaComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/

	switch(this->law){
		case 1:
			GetAlphaViscousComplement(palpha_complement,gauss);
			break;
		case 3:
			GetAlphaHydroComplement(palpha_complement,gauss);
			break;
		case 4:
			GetAlphaTempComplement(palpha_complement,gauss);
			break;
	  default:
			_error_("not supported");
	}

}/*}}}*/
void Friction::GetAlphaHydroComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/

	/*diverse: */
	int         CoupledFlag;
	IssmDouble  q_exp;
	IssmDouble  C_param;
	IssmDouble  As;
	IssmDouble  Neff;
	IssmDouble  n;
	IssmDouble  alpha;
	IssmDouble  Chi,Gamma;
	IssmDouble  vx,vy,vz,vmag;
	IssmDouble  alpha_complement;

	/*Recover parameters: */
	element->GetInputValue(&q_exp,FrictionQEnum);
	element->GetInputValue(&C_param,FrictionCEnum);

	element->GetInputValue(&As,gauss,FrictionAsEnum);
	element->GetInputValue(&n,gauss,MaterialsRheologyNEnum);
	element->parameters->FindParam(&CoupledFlag,FrictionCouplingEnum);

	if (CoupledFlag==1){
		element->GetInputValue(&Neff,gauss,EffectivePressureEnum);
	}
	else{
		element->GetInputValue(&Neff,gauss,FrictionEffectivePressureEnum);
	}

	if(Neff<0)Neff=0;

	//We need the velocity magnitude to evaluate the basal stress:
	switch(dim){
		case 1:
			element->GetInputValue(&vx,gauss,VxEnum);
			vmag=sqrt(vx*vx);
			break;
		case 2:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			vmag=sqrt(vx*vx+vy*vy);
			break;
		case 3:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			element->GetInputValue(&vz,gauss,VzEnum);
			vmag=sqrt(vx*vx+vy*vy+vz*vz);
			break;
		default:
			_error_("not supported");
	}
	//	vmag=100./(3600.*24.*365.);

	if (q_exp==1){
		alpha=1;
	}
	else{
		alpha=(pow(q_exp-1,q_exp-1))/pow(q_exp,q_exp);
	}
	Chi   = vmag/(pow(C_param,n)*pow(Neff,n)*As);
	Gamma = (Chi/(1.+alpha*pow(Chi,q_exp)));
	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0.) alpha_complement=0.;
	else	if(Neff==0.) alpha_complement=0.;
	else	alpha_complement=-(C_param*Neff/(n*vmag)) *
					pow(Gamma,((1.-n)/n)) *
					(Gamma/As - (alpha*q_exp*pow(Chi,q_exp-1.)* Gamma * Gamma/As));

	_assert_(!xIsNan<IssmDouble>(alpha_complement));
	/*Assign output pointers:*/
	*palpha_complement=alpha_complement;
}
/*}}}*/
void Friction::GetAlphaTempComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/
	/*Here, we want to parameterize the friction as a function of temperature
	 *
	 * alpha2 = alpha2_viscous * 1/f(T)
	 *
	 * where f(T) = exp((T-Tpmp)/gamma)
	 */

	/*Intermediaries: */
	IssmDouble  f,T,pressure,Tpmp,gamma;
	IssmDouble  alpha_complement;

	/*Get viscous part*/
	this->GetAlphaViscousComplement(&alpha_complement,gauss);

	/*Get pressure melting point (Tpmp) for local pressure and get current temperature*/
	element->GetInputValue(&T,gauss,TemperatureEnum);
	element->GetInputValue(&pressure,gauss,PressureEnum);
	Tpmp = element->TMeltingPoint(pressure);

	/*Compute scaling parameter*/
	element->parameters->FindParam(&gamma,FrictionGammaEnum);
	alpha_complement = alpha_complement/ exp((T-Tpmp)/gamma);

	/*Assign output pointers:*/
	*palpha_complement=alpha_complement;
}/*}}}*/
void Friction::GetAlphaViscousComplement(IssmDouble* palpha_complement, Gauss* gauss){/*{{{*/

	/* FrictionGetAlpha2 computes alpha2= drag^2 * Neff ^r * vel ^s, with Neff=rho_ice*g*thickness+rho_ice*g*bed, r=q/p and s=1/p. 
	 * FrictionGetAlphaComplement is used in control methods on drag, and it computes: 
	 * alpha_complement= Neff ^r * vel ^s*/

	/*diverse: */
	IssmDouble  r,s;
	IssmDouble  vx,vy,vz,vmag;
	IssmDouble  drag_p,drag_q;
	IssmDouble  Neff;
	IssmDouble  drag_coefficient;
	IssmDouble  bed,thickness,sealevel;
	IssmDouble  alpha_complement;

	/*Recover parameters: */
	element->GetInputValue(&drag_p,FrictionPEnum);
	element->GetInputValue(&drag_q,FrictionQEnum);
	element->GetInputValue(&thickness, gauss,ThicknessEnum);
	element->GetInputValue(&bed, gauss,BaseEnum);
	element->GetInputValue(&sealevel, gauss,SealevelEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);
	IssmDouble rho_water   = element->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice     = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble gravity     = element->GetMaterialParameter(ConstantsGEnum);

	//compute r and q coefficients: */
	r=drag_q/drag_p;
	s=1./drag_p;

	//From bed and thickness, compute effective pressure when drag is viscous:
	Neff=gravity*(rho_ice*thickness+rho_water*(bed-sealevel));
	if(Neff<0)Neff=0;

	//We need the velocity magnitude to evaluate the basal stress:
	switch(dim){
		case 1:
			element->GetInputValue(&vx,gauss,VxEnum);
			vmag=sqrt(vx*vx);
			break;
		case 2:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			vmag=sqrt(vx*vx+vy*vy);
			break;
		case 3:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			element->GetInputValue(&vz,gauss,VzEnum);
			vmag=sqrt(vx*vx+vy*vy+vz*vz);
			break;
		default:
			_error_("not supported");
	}

	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0. && (s-1.)<0.) alpha_complement=0.;
	else alpha_complement=pow(Neff,r)*pow(vmag,(s-1));_assert_(!xIsNan<IssmDouble>(alpha_complement));

	/*Assign output pointers:*/
	*palpha_complement=alpha_complement;
}
/*}}}*/
void Friction::GetAlpha2(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	switch(this->law){
		case 1:
			GetAlpha2Viscous(palpha2,gauss);
			break;
		case 2:
			GetAlpha2Weertman(palpha2,gauss);
			break;
		case 3:
			GetAlpha2Hydro(palpha2,gauss);
			break;
		case 4:
			GetAlpha2Temp(palpha2,gauss);
			break;
		case 5:
			GetAlpha2WaterLayer(palpha2,gauss);
			break;
		case 6:
			GetAlpha2WeertmanTemp(palpha2,gauss);
			break;
		case 7:
			GetAlpha2Coulomb(palpha2,gauss);
			break;
		case 8:
			GetAlpha2Sommers(palpha2,gauss);
			break;
		case 9:
			GetAlpha2Josh(palpha2,gauss);
			break;
	  default:
			_error_("Friction law "<< this->law <<" not supported");
	}

}/*}}}*/
void Friction::GetAlpha2Coulomb(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient 
	  alpha2= drag^2 * Neff ^r * | vel | ^(s-1), with Neff=rho_ice*g*thickness+rho_ice*g*base, r=q/p and s=1/p**/

	/*diverse: */
	IssmDouble  r,s;
	IssmDouble  drag_p, drag_q;
	IssmDouble  Neff;
	IssmDouble  thickness,base,bed,floatation_thickness,sealevel;
	IssmDouble  vx,vy,vz,vmag;
	IssmDouble  drag_coefficient,drag_coefficient_coulomb;
	IssmDouble  alpha2,alpha2_coulomb;

	/*Recover parameters: */
	element->GetInputValue(&drag_p,FrictionPEnum);
	element->GetInputValue(&drag_q,FrictionQEnum);
	element->GetInputValue(&thickness, gauss,ThicknessEnum);
	element->GetInputValue(&base, gauss,BaseEnum);
	element->GetInputValue(&sealevel, gauss,SealevelEnum);
	element->GetInputValue(&bed, gauss,BedEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);
	element->GetInputValue(&drag_coefficient_coulomb, gauss,FrictionCoefficientcoulombEnum);
	IssmDouble rho_water        = element->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice          = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble gravity          = element->GetMaterialParameter(ConstantsGEnum);

	//compute r and q coefficients: */
	r=drag_q/drag_p;
	s=1./drag_p;

	//From base and thickness, compute effective pressure when drag is viscous:
	Neff=gravity*(rho_ice*thickness+rho_water*(base-sealevel));
	if(Neff<0)Neff=0;

	switch(dim){
		case 1:
			element->GetInputValue(&vx,gauss,VxEnum);
			vmag=sqrt(vx*vx);
			break;
		case 2:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			vmag=sqrt(vx*vx+vy*vy);
			break;
		case 3:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			element->GetInputValue(&vz,gauss,VzEnum);
			vmag=sqrt(vx*vx+vy*vy+vz*vz);
			break;
		default:
			_error_("not supported");
	}

	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0. && (s-1.)<0.) alpha2=0.;
	else alpha2=drag_coefficient*drag_coefficient*pow(Neff,r)*pow(vmag,(s-1.));

	floatation_thickness=0;
	if(bed<0) floatation_thickness=-rho_water/rho_ice*bed;
	if(vmag==0.) alpha2_coulomb=0.;
	else alpha2_coulomb=drag_coefficient_coulomb*drag_coefficient_coulomb*rho_water*gravity*(thickness-floatation_thickness)/vmag;

	if(alpha2_coulomb<alpha2) alpha2=alpha2_coulomb;

	_assert_(!xIsNan<IssmDouble>(alpha2));

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Hydro(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient 
		Based on Gagliardini 2007, needs a good effective pressure computation
		Not tested so far so use at your own risks
	  alpha2= NeffC[Chi/(1+alpha*Chi^q)]^(1/n)*1/vel  with
		-Chi=|vel|/(C^n*Neff^n*As)
		-alpha=(q-1)^(q-1)/q^q */

	/*diverse: */
	int         CoupledFlag;
	IssmDouble  q_exp;
	IssmDouble  C_param;
	IssmDouble  As;

	IssmDouble  Neff;
	IssmDouble  n;

	IssmDouble  alpha;
	IssmDouble  Chi,Gamma;

	IssmDouble  vx,vy,vz,vmag;
	IssmDouble  alpha2;

	/*Recover parameters: */
	element->GetInputValue(&q_exp,FrictionQEnum);
	element->GetInputValue(&C_param,FrictionCEnum);
	element->GetInputValue(&As,gauss,FrictionAsEnum);
	element->GetInputValue(&n,gauss,MaterialsRheologyNEnum);
	
	element->parameters->FindParam(&CoupledFlag,FrictionCouplingEnum);
	if (CoupledFlag==1){
		element->GetInputValue(&Neff,gauss,EffectivePressureEnum);
	}
	else{
		element->GetInputValue(&Neff,gauss,FrictionEffectivePressureEnum);
	}
		
	if(Neff<0)Neff=0;

	switch(dim){
		case 1:
			element->GetInputValue(&vx,gauss,VxEnum);
			vmag=sqrt(vx*vx);
			break;
		case 2:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			vmag=sqrt(vx*vx+vy*vy);
			break;
		case 3:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			element->GetInputValue(&vz,gauss,VzEnum);
			vmag=sqrt(vx*vx+vy*vy+vz*vz);
			break;
		default:
			_error_("not supported");
	}

	//	vmag=100./(3600.*24.*365.);
	//compute alpha and Chi coefficients: */
	if (q_exp==1){
		alpha=1;
	}
	else{
		alpha=(pow(q_exp-1,q_exp-1))/pow(q_exp,q_exp);
	}
	Chi=vmag/(pow(C_param,n)*pow(Neff,n)*As);
	Gamma=(Chi/(1. + alpha * pow(Chi,q_exp)));
	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0.) alpha2=0.; 
	else	if (Neff==0) alpha2=0.0;
	else	alpha2=Neff * C_param * pow(Gamma,1./n) * 1/vmag;

	_assert_(!xIsNan<IssmDouble>(alpha2));
	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Sommers(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/* FrictionGetAlpha2 computes alpha2= drag^2 * Neff, with Neff=rho_ice*g*thickness+rho_ice*g*(head-bed)*/

	/*diverse: */
	IssmDouble  pressure_ice,pressure_water;
	IssmDouble  Neff;
	IssmDouble  drag_coefficient;
	IssmDouble  bed,thickness,head,sealevel;
	IssmDouble  alpha2;

	/*Recover parameters: */
	element->GetInputValue(&thickness, gauss,ThicknessEnum);
	element->GetInputValue(&bed, gauss,BaseEnum);
	element->GetInputValue(&head, gauss,HydrologyHeadEnum);
	element->GetInputValue(&sealevel, gauss,SealevelEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);
	IssmDouble rho_water   = element->GetMaterialParameter(MaterialsRhoFreshwaterEnum);
	IssmDouble rho_ice     = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble gravity     = element->GetMaterialParameter(ConstantsGEnum);

	//From bed and thickness, compute effective pressure when drag is viscous:
	pressure_ice   = rho_ice*gravity*thickness;
	pressure_water = rho_water*gravity*(head-bed+sealevel);
	Neff=pressure_ice-pressure_water;
	if(Neff<0.) Neff=0.;

	alpha2=drag_coefficient*drag_coefficient*Neff;
	_assert_(!xIsNan<IssmDouble>(alpha2));

	/*Assign output pointers:*/
	*palpha2=alpha2;
}
/*}}}*/
void Friction::GetAlpha2Temp(IssmDouble* palpha2, Gauss* gauss){/*{{{*/
	/*Here, we want to parameterize the friction as a function of temperature
	 *
	 * alpha2 = alpha2_viscous * 1/f(T)
	 *
	 * where f(T) = exp((T-Tpmp)/gamma)
	 */

	/*Intermediaries: */
	IssmDouble  f,T,pressure,Tpmp,gamma;
	IssmDouble  alpha2;

	/*Get viscous part*/
	this->GetAlpha2Viscous(&alpha2,gauss);

	/*Get pressure melting point (Tpmp) for local pressure and get current temperature*/
	element->GetInputValue(&T,gauss,TemperatureEnum);
	element->GetInputValue(&pressure,gauss,PressureEnum);
	Tpmp = element->TMeltingPoint(pressure);

	/*Compute scaling parameter*/
	element->parameters->FindParam(&gamma,FrictionGammaEnum);
	alpha2 = alpha2 / exp((T-Tpmp)/gamma);

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Josh(IssmDouble* palpha2, Gauss* gauss){/*{{{*/
	/*Here, we want to parameterize the friction as a function of temperature
	 *
	 * alpha2 = alpha2_viscous * 1/f(T)
	 *
	 * where f(T) = exp((T-Tpmp)/gamma)
	 */

	/*Intermediaries: */
	IssmDouble  T,Tpmp,deltaT,deltaTref,pressure;
	IssmDouble  alpha2,time,gamma;
	const IssmDouble yts = 365*24*3600.;

	/*Get viscous part*/
	this->GetAlpha2Viscous(&alpha2,gauss);

	/*Get delta Refs*/
	element->GetInputValue(&deltaTref,gauss,FrictionPressureAdjustedTemperatureEnum);

	/*Compute delta T*/
	element->GetInputValue(&T,gauss,TemperatureEnum);
	element->GetInputValue(&pressure,gauss,PressureEnum);
	Tpmp = element->TMeltingPoint(pressure);
	deltaT = T-Tpmp;

	/*Compute gamma*/
	element->parameters->FindParam(&time,TimeEnum);
	element->parameters->FindParam(&gamma,FrictionGammaEnum);
	//if(time<25e3*yts){
	//	gamma = 10.;
	//}
	//else{
	//	gamma = 5.;
	//}
	//gamma = 5.;

	/*Compute scaling parameter*/
	alpha2 = alpha2 * exp((deltaTref - deltaT)/(2*gamma));

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Viscous(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient 
	  alpha2= drag^2 * Neff ^r * | vel | ^(s-1), with Neff=rho_ice*g*thickness+rho_ice*g*bed, r=q/p and s=1/p**/

	/*diverse: */
	IssmDouble  r,s;
	IssmDouble  drag_p, drag_q;
	IssmDouble  Neff;
	IssmDouble  thickness,base,sealevel;
	IssmDouble  vx,vy,vz,vmag;
	IssmDouble  drag_coefficient;
	IssmDouble  alpha2;

	/*Recover parameters: */
	element->GetInputValue(&drag_p,FrictionPEnum);
	element->GetInputValue(&drag_q,FrictionQEnum);
	element->GetInputValue(&thickness, gauss,ThicknessEnum);
	element->GetInputValue(&base, gauss,BaseEnum);
	element->GetInputValue(&sealevel, gauss,SealevelEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);
	IssmDouble rho_water   = element->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice     = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble gravity     = element->GetMaterialParameter(ConstantsGEnum);

	//compute r and q coefficients: */
	r=drag_q/drag_p;
	s=1./drag_p;

	//From base and thickness, compute effective pressure when drag is viscous:
	Neff=gravity*(rho_ice*thickness+rho_water*(base-sealevel));
	if(Neff<0)Neff=0;

	switch(dim){
		case 1:
			element->GetInputValue(&vx,gauss,VxEnum);
			vmag=sqrt(vx*vx);
			break;
		case 2:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			vmag=sqrt(vx*vx+vy*vy);
			break;
		case 3:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			element->GetInputValue(&vz,gauss,VzEnum);
			vmag=sqrt(vx*vx+vy*vy+vz*vz);
			break;
		default:
			_error_("not supported");
	}

	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0. && (s-1.)<0.) alpha2=0.;
	else alpha2=drag_coefficient*drag_coefficient*pow(Neff,r)*pow(vmag,(s-1.));
	_assert_(!xIsNan<IssmDouble>(alpha2));

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2WaterLayer(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient 
	  alpha2= drag^2 * Neff ^r * | vel | ^(s-1), with Neff=rho_ice*g*thickness+rho_ice*g*bed, r=q/p and s=1/p**/

	/*diverse: */
	IssmDouble  r,s;
	IssmDouble  drag_p, drag_q;
	IssmDouble  Neff,F;
	IssmDouble  thickness,bed,sealevel;
	IssmDouble  vx,vy,vz,vmag;
	IssmDouble  drag_coefficient,water_layer;
	IssmDouble  alpha2;

	/*Recover parameters: */
	element->parameters->FindParam(&F,FrictionFEnum);
	element->GetInputValue(&drag_p,FrictionPEnum);
	element->GetInputValue(&drag_q,FrictionQEnum);
	element->GetInputValue(&thickness, gauss,ThicknessEnum);
	element->GetInputValue(&bed, gauss,BaseEnum);
	element->GetInputValue(&sealevel, gauss,SealevelEnum);
	element->GetInputValue(&drag_coefficient, gauss,FrictionCoefficientEnum);
	element->GetInputValue(&water_layer, gauss,FrictionWaterLayerEnum);
	IssmDouble rho_water   = element->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice     = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble gravity     = element->GetMaterialParameter(ConstantsGEnum);

	//compute r and q coefficients: */
	r=drag_q/drag_p;
	s=1./drag_p;

	//From bed and thickness, compute effective pressure when drag is viscous:
	if(bed>0) bed=0;
	if(water_layer==0) Neff=gravity*rho_ice*thickness+gravity*rho_water*(bed-sealevel);
	else if(water_layer>0) Neff=gravity*rho_ice*thickness*F;
	else _error_("negative water layer thickness");
	if(Neff<0) Neff=0;

	switch(dim){
		case 1:
			element->GetInputValue(&vx,gauss,VxEnum);
			vmag=sqrt(vx*vx);
			break;
		case 2:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			vmag=sqrt(vx*vx+vy*vy);
			break;
		case 3:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			element->GetInputValue(&vz,gauss,VzEnum);
			vmag=sqrt(vx*vx+vy*vy+vz*vz);
			break;
		default:
			_error_("not supported");
	}

	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0. && (s-1.)<0.) alpha2=0.;
	else alpha2=drag_coefficient*drag_coefficient*pow(Neff,r)*pow(vmag,(s-1.));
	_assert_(!xIsNan<IssmDouble>(alpha2));

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2Weertman(IssmDouble* palpha2, Gauss* gauss){/*{{{*/

	/*This routine calculates the basal friction coefficient alpha2= C^-1/m |v|^(1/m-1) */

	/*diverse: */
	IssmDouble  C,m;
	IssmDouble  vx,vy,vz,vmag;
	IssmDouble  alpha2;

	/*Recover parameters: */
	element->GetInputValue(&C,gauss,FrictionCEnum);
	element->GetInputValue(&m,FrictionMEnum);

	switch(dim){
		case 1:
			element->GetInputValue(&vx,gauss,VxEnum);
			vmag=sqrt(vx*vx);
			break;
		case 2:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			vmag=sqrt(vx*vx+vy*vy);
			break;
		case 3:
			element->GetInputValue(&vx,gauss,VxEnum);
			element->GetInputValue(&vy,gauss,VyEnum);
			element->GetInputValue(&vz,gauss,VzEnum);
			vmag=sqrt(vx*vx+vy*vy+vz*vz);
			break;
		default:
			_error_("not supported");
	}

	/*Check to prevent dividing by zero if vmag==0*/
	if(vmag==0. && (1./m-1.)<0.) alpha2=0.;
	else alpha2=pow(C,-1./m)*pow(vmag,(1./m-1.));
	_assert_(!xIsNan<IssmDouble>(alpha2));

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
void Friction::GetAlpha2WeertmanTemp(IssmDouble* palpha2, Gauss* gauss){/*{{{*/
	/*Here, we want to parameterize the friction as a function of temperature
	 *
	 * alpha2 = alpha2_weertman * 1/f(T)
	 *
	 * where f(T) = exp((T-Tpmp)/gamma)
	 */

	/*Intermediaries: */
	IssmDouble  f,T,pressure,Tpmp,gamma;
	IssmDouble  alpha2;

	/*Get viscous part*/
	this->GetAlpha2Weertman(&alpha2,gauss);

	/*Get pressure melting point (Tpmp) for local pressure and get current temperature*/
	element->GetInputValue(&T,gauss,TemperatureEnum);
	element->GetInputValue(&pressure,gauss,PressureEnum);
	Tpmp = element->TMeltingPoint(pressure);

	/*Compute scaling parameter*/
	element->parameters->FindParam(&gamma,FrictionGammaEnum);
	alpha2 = alpha2 / exp((T-Tpmp)/gamma);

	/*Assign output pointers:*/
	*palpha2=alpha2;
}/*}}}*/
