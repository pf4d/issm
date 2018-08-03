/*!\file Matice.c
 * \brief: implementation of the Matice object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Matice.h"
#include "./Materials.h"
#include "../Inputs/Input.h"
#include "../Inputs/Inputs.h"
#include "../Inputs/TriaInput.h"
#include "../Inputs/PentaInput.h"
#include "../Inputs/ControlInput.h"
#include "../Elements/Element.h"
#include "../Elements/Tria.h"
#include "../Elements/Penta.h"
#include "../Params/Parameters.h"
#include "../Vertex.h"
#include "../Hook.h"
#include "../Node.h"
#include "../IoModel.h"
#include "../../shared/shared.h"

/*Matice constructors and destructor*/
Matice::Matice(){/*{{{*/
	this->helement=NULL;
	this->element=NULL;
	this->isdamaged=false;
	this->isenhanced=false;
	return;
}
/*}}}*/
Matice::Matice(int matice_mid,int index, IoModel* iomodel){/*{{{*/

	 /*Get material type and initialize object*/
   int materialtype;
   iomodel->FindConstant(&materialtype,"md.materials.type");
	this->Init(matice_mid,index,materialtype);

}
/*}}}*/
Matice::Matice(int matice_mid,int index,int materialtype){/*{{{*/

	this->Init(matice_mid,index,materialtype);
	return;
} /*}}}*/
Matice::~Matice(){/*{{{*/
	delete helement;
	return;
}
/*}}}*/
void Matice::Init(int matice_mid,int index,int materialtype){/*{{{*/

	/*Initialize id*/
	this->mid=matice_mid;

	/*Hooks: */
	int matice_eid=index+1;
	this->helement=new Hook(&matice_eid,1);
	this->element=NULL;

	/*Material specific properties*/
	switch(materialtype){
		case MatdamageiceEnum:
			this->isdamaged = true;
			this->isenhanced = false;
			break;
		case MaticeEnum:
			this->isdamaged = false;
			this->isenhanced = false;
			break;
		case MatenhancediceEnum:
			this->isdamaged = false;
			this->isenhanced = true;
			break;
		default:
			_error_("Material type not recognized");
	}

	return;
} /*}}}*/

/*Object virtual functions definitions:*/
Object*   Matice::copy() {/*{{{*/

	/*Output*/
	Matice* matice=NULL;

	/*Initialize output*/
	matice=new Matice();

	/*copy fields: */
	matice->mid=this->mid;
	matice->helement=(Hook*)this->helement->copy();
	matice->element =(Element*)this->helement->delivers();
	matice->isdamaged = this->isdamaged;
	matice->isenhanced = this->isenhanced;

	return matice;
}
/*}}}*/
Material* Matice::copy2(Element* element_in) {/*{{{*/

	/*Output*/
	Matice* matice=NULL;

	/*Initialize output*/
	matice=new Matice();

	/*copy fields: */
	matice->mid=this->mid;
	matice->helement=(Hook*)this->helement->copy();
	matice->element =element_in;
	matice->isdamaged = this->isdamaged;
	matice->isenhanced = this->isenhanced;

	return matice;
}
/*}}}*/
void      Matice::DeepEcho(void){/*{{{*/

	_printf_("Matice:\n");
	_printf_("   mid: " << mid << "\n");
	_printf_("   isdamaged: " << isdamaged << "\n");
	_printf_("   isenhanced: " << isenhanced << "\n");

	/*helement and element DeepEcho were commented to avoid recursion.*/
	/*Example: element->DeepEcho calls matice->DeepEcho which calls element->DeepEcho etc*/
	_printf_("   helement:\n");
	_printf_("		note: helement not printed to avoid recursion.\n");
	//if(helement) helement->DeepEcho();
	//else _printf_("   helement = NULL\n");
	
	_printf_("   element:\n");
	_printf_("     note: element not printed to avoid recursion.\n");
	//if(element) element->DeepEcho();
	//else _printf_("   element = NULL\n");
}		
/*}}}*/
void      Matice::Echo(void){/*{{{*/
	
	_printf_("Matice:\n");
	_printf_("   mid: " << mid << "\n");
	_printf_("   isdamaged: " << isdamaged << "\n");
	_printf_("   isenhanced: " << isenhanced << "\n");
	
	/*helement and element Echo were commented to avoid recursion.*/
	/*Example: element->Echo calls matice->Echo which calls element->Echo etc*/
	_printf_("   helement:\n");
	_printf_("     note: helement not printed to avoid recursion.\n");
	//if(helement) helement->Echo();
	//else _printf_("   helement = NULL\n");
	
	_printf_("   element:\n");
	_printf_("     note: element not printed to avoid recursion.\n");
	//if(element) element->Echo();
	//else _printf_("   element = NULL\n");
}
/*}}}*/
int       Matice::Id(void){ return mid; }/*{{{*/
/*}}}*/
void      Matice::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	if(marshall_direction==MARSHALLING_BACKWARD)helement=new Hook(); 
	
	MARSHALLING_ENUM(MaticeEnum);
	MARSHALLING(mid);
	MARSHALLING(isdamaged);
	MARSHALLING(isenhanced);
	this->helement->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->element=(Element*)this->helement->delivers();

}
/*}}}*/
int       Matice::ObjectEnum(void){/*{{{*/

	return MaticeEnum;

}
/*}}}*/

/*Matice management*/
void  Matice::Configure(Elements* elementsin){/*{{{*/

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective 
	 * datasets, using internal ids and offsets hidden in hooks: */
	helement->configure((DataSet*)elementsin);
	this->element  = (Element*)helement->delivers();
}
/*}}}*/
IssmDouble Matice::GetA(){/*{{{*/
	/*
	 * A = 1/B^n
	 */

	IssmDouble B,n;

	element->inputs->GetInputAverage(&B,MaterialsRheologyBEnum);
	n=this->GetN();

	return pow(B,-n);
}
/*}}}*/
IssmDouble Matice::GetAbar(){/*{{{*/
	/*
	 * A = 1/B^n
	 */

	IssmDouble B,n;

	element->inputs->GetInputAverage(&B,MaterialsRheologyBbarEnum);
	n=this->GetN();

	return pow(B,-n);
}
/*}}}*/
IssmDouble Matice::GetB(){/*{{{*/

	/*Output*/
	IssmDouble B;

	element->inputs->GetInputAverage(&B,MaterialsRheologyBEnum);
	return B;
}
/*}}}*/
IssmDouble Matice::GetBbar(){/*{{{*/

	/*Output*/
	IssmDouble Bbar;

	element->inputs->GetInputAverage(&Bbar,MaterialsRheologyBbarEnum);
	return Bbar;
}
/*}}}*/
IssmDouble Matice::GetD(){/*{{{*/

	_assert_(this->isdamaged);
	/*Output*/
	IssmDouble D;
	if(this->isdamaged)element->inputs->GetInputAverage(&D,DamageDEnum);
	return D;
}
/*}}}*/
IssmDouble Matice::GetDbar(){/*{{{*/

	_assert_(this->isdamaged);
	/*Output*/
	IssmDouble Dbar;
	if(this->isdamaged)element->inputs->GetInputAverage(&Dbar,DamageDbarEnum);
	return Dbar;
}
/*}}}*/
IssmDouble Matice::GetE(){/*{{{*/

	_assert_(this->isenhanced);
	/*Output*/
	IssmDouble E;
	if(this->isenhanced)element->inputs->GetInputAverage(&E,MaterialsRheologyEEnum);
	return E;
}
/*}}}*/
IssmDouble Matice::GetEbar(){/*{{{*/

	_assert_(this->isenhanced);
	/*Output*/
	IssmDouble Ebar;
	if(this->isenhanced)element->inputs->GetInputAverage(&Ebar,MaterialsRheologyEbarEnum);
	return Ebar;
}
/*}}}*/
IssmDouble Matice::GetN(){/*{{{*/

	/*Output*/
	IssmDouble n;

	element->inputs->GetInputAverage(&n,MaterialsRheologyNEnum);
	return n;
}
/*}}}*/
bool Matice::IsDamage(){/*{{{*/

	return this->isdamaged;
}
/*}}}*/
bool Matice::IsEnhanced(){/*{{{*/

	return this->isenhanced;
}
/*}}}*/
void  Matice::GetViscosity(IssmDouble* pviscosity,IssmDouble eps_eff){/*{{{*/
	/*From a string tensor and a material object, return viscosity, using Glen's flow law.
								(1-D) B
	  viscosity= -------------------------
						  2 E^[1/n] eps_eff ^[(n-1)/n]

	  where viscosity is the viscosity, B the flow law parameter , eps_eff is the effective strain rate,
	  n the flow law exponent, and E is the enhancement factor.

	  If eps_eff = 0 , it means this is the first time SystemMatrices is being run, and we 
	  return 10^14, initial viscosity.
	  */

	/*output: */
	IssmDouble viscosity;

	/*Intermediary: */
	IssmDouble B,D=0.,E=1.,n;

	/*Get B and n*/
	B=GetB(); _assert_(B>0.);
	n=GetN(); _assert_(n>0.);
	if(this->isdamaged){
		D=GetD();
		_assert_(D>=0. && D<1.);
	}
	if(this->isenhanced){
		E=GetE();
		_assert_(E>0.);
	}

	if (n==1.){
		/*Linear Viscous behavior (Newtonian fluid) viscosity=B/2E: */
		viscosity=(1.-D)*B/(2.*E);
	}
	else{

		/*if no strain rate, return maximum viscosity*/
		if(eps_eff==0.){
			viscosity = 1.e+14/2.;
			//viscosity = B;
			//viscosity=2.5*pow(10.,17);
		}

		else{
			viscosity=(1.-D)*B/(2.*pow(E,1./n)*pow(eps_eff,(n-1.)/n));
		}
	}

	/*Checks in debugging mode*/
	if(viscosity<=0) _error_("Negative viscosity");

	/*Return: */
	*pviscosity=viscosity;
}
/*}}}*/
void  Matice::GetViscosityBar(IssmDouble* pviscosity,IssmDouble eps_eff){/*{{{*/
	/*From a string tensor and a material object, return viscosity, using Glen's flow law.
								(1-D) B
	  viscosity= -------------------------
						  2 E^[1/n] eps_eff ^[(n-1)/n]

	  where B the flow law parameter, eps_eff is the effective strain rate, n the flow law exponent,
	  and E is the enhancement factor.

	  If eps_eff = 0 , it means this is the first time SystemMatrices is being run, and we 
	  return 10^14, initial viscosity.
	  */

	/*output: */
	IssmDouble viscosity;

	/*Intermediary: */
	IssmDouble B,D=0.,E=1.,n;

	/*Get B and n*/
	B=GetBbar(); _assert_(B>0.);
	n=GetN();    _assert_(n>0.);
	if(this->isdamaged){
		D=GetDbar();
		_assert_(D>=0. && D<1.);
	}
	if(this->isenhanced){
		E=GetEbar();
		_assert_(E>0.);
	}

	if (n==1.){
		/*Linear Viscous behavior (Newtonian fluid) viscosity=B/2E: */
		viscosity=(1.-D)*B/(2.*E);
	}
	else{

		/*if no strain rate, return maximum viscosity*/
		if(eps_eff==0.){
			viscosity = 1.e+14/2.;
			//viscosity=2.5*pow(10.,17);
		}

		else{
			viscosity=(1.-D)*B/(2.*pow(E,1./n)*pow(eps_eff,(n-1.)/n));
		}
	}

	/*Checks in debugging mode*/
	if(viscosity<=0) _error_("Negative viscosity");

	/*Return: */
	*pviscosity=viscosity;
}
/*}}}*/
void  Matice::GetViscosityComplement(IssmDouble* pviscosity_complement, IssmDouble* epsilon){/*{{{*/
	/*Return viscosity accounting for steady state power law creep [Thomas and SSA, 1982]: 
	 *
	 *  										                (1-D)
	 * viscosity= -------------------------------------------------------------------
	 *  				  2[ exx^2+eyy^2+exx*eyy+exy^2+exz^2+eyz^2 ]^[(n-1)/2n]
	 *
	 * If epsilon is NULL, it means this is the first time Gradjb is being run, and we 
	 * return mu20, initial viscosity.
	 */

	/*output: */
	IssmDouble viscosity_complement;

	/*input strain rate: */
	IssmDouble exx,eyy,exy;

	/*Intermediary value A and exponent e: */
	IssmDouble A,e;
	IssmDouble D=0.,n;

	/*Get D and n*/
	if(this->isdamaged){
		D=GetDbar(); /* GetD()? */
		_assert_(D>=0. && D<1.);
	}
	n=GetN();

	if(epsilon){
		exx=*(epsilon+0);
		eyy=*(epsilon+1);
		exy=*(epsilon+2);

		/*Build viscosity: mu2=(1-D)/(2*A^e) */
		A=pow(exx,2)+pow(eyy,2)+pow(exy,2)+exx*eyy;
		if(A==0){
			/*Maximum viscosity_complement for 0 shear areas: */
			viscosity_complement=2.25*pow(10.,17);
		}
		else{
			e=(n-1)/(2*n);

			viscosity_complement=(1-D)/(2*pow(A,e));
		}
	}
	else{
		viscosity_complement=4.5*pow(10.,17);
	}

	/*Checks in debugging mode*/
	_assert_(D>=0 && D<1);
	_assert_(n>0);
	_assert_(viscosity_complement>0);

	/*Return: */
	*pviscosity_complement=viscosity_complement;
}
/*}}}*/
void  Matice::GetViscosityDComplement(IssmDouble* pviscosity_complement, IssmDouble* epsilon){/*{{{*/
	/*Return viscosity derivative for control method d(mu)/dD: 
	 *
	 *  										               B 
	 * dviscosity= - -------------------------------------------------------------------
	 *  				  2[ exx^2+eyy^2+exx*eyy+exy^2+exz^2+eyz^2 ]^[(n-1)/2n]
	 *
	 * If epsilon is NULL, it means this is the first time Gradjb is being run, and we 
	 * return mu20, initial viscosity.
	 */

	/*output: */
	IssmDouble viscosity_complement;

	/*input strain rate: */
	IssmDouble exx,eyy,exy;

	/*Intermediary value A and exponent e: */
	IssmDouble A,e;
	IssmDouble B,n;

	/*Get B and n*/
	B=GetBbar();
	n=GetN();

	if(epsilon){
		exx=*(epsilon+0);
		eyy=*(epsilon+1);
		exy=*(epsilon+2);

		/*Build viscosity: mu2=B/(2*A^e) */
		A=pow(exx,2)+pow(eyy,2)+pow(exy,2)+exx*eyy;
		if(A==0){
			/*Maximum viscosity_complement for 0 shear areas: */
			viscosity_complement=- 2.25*pow(10.,17);
		}
		else{
			e=(n-1)/(2*n);

			viscosity_complement=- B/(2*pow(A,e));
		}
	}
	else{
		viscosity_complement=- 4.5*pow(10.,17);
	}

	/*Checks in debugging mode*/
	_assert_(B>0);
	_assert_(n>0);
	_assert_(viscosity_complement<0);

	/*Return: */
	*pviscosity_complement=viscosity_complement;
}
/*}}}*/
void  Matice::GetViscosityDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* epsilon){/*{{{*/

	/*output: */
	IssmDouble mu_prime;
	IssmDouble mu,n,eff2;

	/*input strain rate: */
	IssmDouble exx,eyy,exy,exz,eyz;


	if((epsilon[0]==0) && (epsilon[1]==0) && (epsilon[2]==0) && 
				(epsilon[3]==0) && (epsilon[4]==0)){
		mu_prime=0.5*pow(10.,14);
	}
	else{

		/*Retrive strain rate components: */
		exx=epsilon[0];
		eyy=epsilon[1];
		exy=epsilon[2];
		exz=epsilon[3];
		eyz=epsilon[4];
		eff2 = exx*exx + eyy*eyy + exx*eyy + exy*exy + exz*exz + eyz*eyz;

		GetViscosity(&mu,sqrt(eff2));
		n=GetN();
		mu_prime=(1.-n)/(2.*n) * mu/eff2;
	}

	/*Assign output pointers:*/
	*pmu_prime=mu_prime;
}
/*}}}*/
void  Matice::GetViscosity_B(IssmDouble* pdmudB,IssmDouble eps_eff){/*{{{*/

	/*output: */
	IssmDouble dmudB;

	/*Intermediary: */
	IssmDouble D=0.,E=1.,n;

	/*Get B and n*/
	n=GetN(); _assert_(n>0.);
	if(this->isdamaged){
		D=GetD();
		_assert_(D>=0. && D<1.);
	}
	if(this->isenhanced){
		E=GetE();
		_assert_(E>0.);
	}

	if(n==1.){
		/*Linear Viscous behavior (Newtonian fluid) dmudB=B/2E: */
		dmudB=(1.-D)/(2.*E);
	}
	else{
		if(eps_eff==0.) dmudB = 0.;
		else            dmudB = (1.-D)/(2.*pow(E,1./n)*pow(eps_eff,(n-1.)/n));
	}

	/*Return: */
	*pdmudB=dmudB;
}
/*}}}*/
void  Matice::GetViscosity_D(IssmDouble* pdmudD,IssmDouble eps_eff){/*{{{*/

	/*output: */
	IssmDouble dmudD;

	/*Intermediary: */
	IssmDouble n,B,E=1.;

	/*Get B and n*/
	n=GetN(); _assert_(n>0.);
	B=GetBbar();
	_assert_(this->isdamaged);
	if(this->isenhanced){
		E=GetE();
		_assert_(E>0.);
	}

	if(n==1.){
		/*Linear Viscous behavior (Newtonian fluid) dmudB=B/2E: */
		dmudD=-B/(2.*E);
	}
	else{
		if(eps_eff==0.) dmudD = 0.;
		else            dmudD = -B/(2.*pow(E,1./n)*pow(eps_eff,(n-1.)/n));
	}

	/*Return: */
	*pdmudD=dmudD;
}
/*}}}*/
void  Matice::GetViscosity2dDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* epsilon){/*{{{*/

	/*output: */
	IssmDouble mu_prime;
	IssmDouble mu,n,eff2;

	/*input strain rate: */
	IssmDouble exx,eyy,exy;

	if((epsilon[0]==0) && (epsilon[1]==0) && (epsilon[2]==0)){
		mu_prime=0.5*pow(10.,14);
	}
	else{
		/*Retrive strain rate components: */
		exx=epsilon[0];
		eyy=epsilon[1];
		exy=epsilon[2];
		eff2 = exx*exx + eyy*eyy + exx*eyy + exy*exy ;

		GetViscosityBar(&mu,sqrt(eff2));
		n=GetN();
		mu_prime=(1.-n)/(2.*n)*mu/eff2;
	}

	/*Assign output pointers:*/
	*pmu_prime=mu_prime;
}
/*}}}*/
void  Matice::InputUpdateFromConstant(IssmDouble constant, int name){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Matice::InputUpdateFromConstant(int constant, int name){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Matice::InputUpdateFromConstant(bool constant, int name){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Matice::InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrows, int ncols,int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Matice::InputUpdateFromVector(IssmDouble* vector, int name, int type){/*{{{*/

}
/*}}}*/
void  Matice::InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){/*{{{*/

}
/*}}}*/
void  Matice::ResetHooks(){/*{{{*/

	this->element=NULL;

	/*Get Element type*/
	this->helement->reset();

}
/*}}}*/
void  Matice::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Matice::ViscosityFS(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/
	/*The effective strain rate is defined in Paterson 3d Ed p 91 eq 9,
	 * and Cuffey p 303 eq 8.18:
	 *
	 *  2 eps_eff^2 = eps_xx^2 + eps_yy^2 + eps_zz^2 + 2(eps_xy^2 + eps_xz^2 + eps_yz^2)
	 *
	 *  or
	 *
	 *  eps_eff = 1/sqrt(2) sqrt( \sum_ij eps_ij^2 )
	 *
	 *          = 1/sqrt(2) ||eps||_F
	 *
	 *  where ||.||_F is the Frobenius norm */

	/*Intermediaries*/
	IssmDouble viscosity;
	IssmDouble epsilon3d[6]; /* epsilon=[exx,eyy,ezz,exy,exz,eyz];*/
	IssmDouble epsilon2d[3]; /* epsilon=[exx,eyy,exy];            */
	IssmDouble eps_eff;
	IssmDouble eps0=1.e-27;

	if(dim==3){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		element->StrainRateFS(&epsilon3d[0],xyz_list,gauss,vx_input,vy_input,vz_input);
		eps_eff = sqrt(epsilon3d[0]*epsilon3d[0] + epsilon3d[1]*epsilon3d[1] + epsilon3d[3]*epsilon3d[3] +  epsilon3d[4]*epsilon3d[4] + epsilon3d[5]*epsilon3d[5] + epsilon3d[0]*epsilon3d[1]+eps0*eps0);
	}
	else{
		/* eps_eff^2 = 1/2 ( exx^2 + eyy^2 + 2*exy^2 )*/
		element->StrainRateSSA(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = 1./sqrt(2.)*sqrt(epsilon2d[0]*epsilon2d[0] + epsilon2d[1]*epsilon2d[1] + 2.*epsilon2d[2]*epsilon2d[2]);
	}

	/*Get viscosity*/
	this->GetViscosity(&viscosity,eps_eff);

	/*Assign output pointer*/
	*pviscosity=viscosity;
}
/*}}}*/
void  Matice::ViscosityFSDerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon){/*{{{*/
	this->GetViscosityDerivativeEpsSquare(pmu_prime,epsilon);
}/*}}}*/
void  Matice::ViscosityHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble viscosity;
	IssmDouble epsilon3d[5];/* epsilon=[exx,eyy,exy,exz,eyz];*/
	IssmDouble epsilon2d[2];/* epsilon=[exx,exy];            */
	IssmDouble eps_eff,E=1.0;

	if(dim==3){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exz^2 + eyz^2 + exx*eyy */
		element->StrainRateHO(&epsilon3d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = sqrt(epsilon3d[0]*epsilon3d[0] + epsilon3d[1]*epsilon3d[1] + epsilon3d[2]*epsilon3d[2] +  epsilon3d[3]*epsilon3d[3] + epsilon3d[4]*epsilon3d[4] + epsilon3d[0]*epsilon3d[1]);
	}
	else{
		/* eps_eff^2 = 1/2 (2*exx^2 + 2*exy^2 ) (since eps_zz = - eps_xx)*/
		element->StrainRateHO2dvertical(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = 1./sqrt(2.)*sqrt(2*epsilon2d[0]*epsilon2d[0] + 2*epsilon2d[1]*epsilon2d[1]);
	}

	/*Get viscosity*/
	this->GetViscosity(&viscosity,eps_eff);
	_assert_(!xIsNan<IssmDouble>(viscosity));

	/*Assign output pointer*/
	*pviscosity=viscosity;
}/*}}}*/
void  Matice::ViscosityHODerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon){/*{{{*/
	this->GetViscosityDerivativeEpsSquare(pmu_prime,epsilon);
}/*}}}*/
void  Matice::ViscosityL1L2(IssmDouble* pviscosity,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* surface_input){/*{{{*/
	/*Compute the L1L2 viscosity
	 *
	 *      1
	 * mu = - A^-1 (sigma'_e)^(1-n)
	 *      2
	 *
	 * sigma'_e^2 = |sigma'_//|^2 + |sigma'_perp|^2 (see Perego 2012 eq. 17,18)
	 *
	 * L1L2 assumptions:
	 *
	 * (1) |eps_b|_// = A (|sigma'_//|^2 + |sigma'_perp|^2)^((n-1)/2) |sigma'_//|
	 * (2) |sigma'_perp|^2 = |rho g (s-z) grad(s)|^2
	 *
	 * Assuming that n = 3, we have a polynom of degree 3 to solve (the only unkown is X=|sigma'_//|)
	 *
	 * A X^3 + A |rho g (s-z) grad(s)|^2 X - |eps_b|_// = 0     */

	IssmDouble z,s,viscosity,p,q,delta;
	IssmDouble tau_perp,tau_par,eps_b,A;
	IssmDouble epsilon[5];   /*exx eyy exy exz eyz*/
	IssmDouble slope[3];

	/*Check that both inputs have been found*/
	if (!vx_input || !vy_input || !surface_input) _error_("Input missing");

	/*Get tau_perp*/
	surface_input->GetInputValue(&s,gauss);
	surface_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);
	z=this->element->GetZcoord(xyz_list,gauss);
	tau_perp = element->matpar->GetMaterialParameter(MaterialsRhoIceEnum) * element->matpar->GetMaterialParameter(ConstantsGEnum) * fabs(s-z)*sqrt(slope[0]*slope[0]+slope[1]*slope[1]);

	/* Get eps_b*/
	element->StrainRateHO(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
	eps_b = sqrt(epsilon[0]*epsilon[0] + epsilon[1]*epsilon[1] + epsilon[0]*epsilon[1] + epsilon[2]*epsilon[2]);
	if(eps_b==0.){
		*pviscosity = 2.5e+17;
		return;
	}

	/*Get A*/
	_assert_(this->GetN()==3.0);
	A=this->GetA();

	/*Solve for tau_perp (http://fr.wikipedia.org/wiki/Méthode_de_Cardan)*/
	p     = tau_perp *tau_perp;
	q     = - eps_b/A;
	delta = q *q + p*p*p*4./27.;
	_assert_(delta>0);
	tau_par = pow(0.5*(-q+sqrt(delta)),1./3.) - pow(0.5*(q+sqrt(delta)),1./3.);

	/*Viscosity*/
	viscosity = 1./(2.*A*(tau_par*tau_par + tau_perp*tau_perp));
	_assert_(!xIsNan(viscosity));
	_assert_(viscosity > 0.);

	/*Assign output pointer*/
	*pviscosity = viscosity;
}/*}}}*/
void  Matice::ViscositySSA(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble viscosity;
	IssmDouble epsilon2d[3];/* epsilon=[exx,eyy,exy];    */
	IssmDouble epsilon1d;   /* epsilon=[exx];    */
	IssmDouble eps_eff;

	if(dim==2){
		/* eps_eff^2 = exx^2 + eyy^2 + exy^2 + exx*eyy*/
		element->StrainRateSSA(&epsilon2d[0],xyz_list,gauss,vx_input,vy_input);
		eps_eff = sqrt(epsilon2d[0]*epsilon2d[0] + epsilon2d[1]*epsilon2d[1] + epsilon2d[2]*epsilon2d[2] + epsilon2d[0]*epsilon2d[1]);
	}
	else{
		/* eps_eff^2 = exx^2*/
		element->StrainRateSSA1d(&epsilon1d,xyz_list,gauss,vx_input);
		eps_eff = fabs(epsilon1d);
	}

	/*Get viscosity*/
	this->GetViscosityBar(&viscosity,eps_eff);

	/*Assign output pointer*/
	*pviscosity=viscosity;
}/*}}}*/
void  Matice::ViscositySSADerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon){/*{{{*/
	this->GetViscosity2dDerivativeEpsSquare(pmu_prime,epsilon);
}/*}}}*/
