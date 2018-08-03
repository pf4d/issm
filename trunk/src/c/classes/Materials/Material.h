/*!\file:  Material.h
 * \brief abstract class for Material object
 */ 

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

/*Headers:*/
/*{{{*/
class Inputs;
template <class doubletype> class Vector;
#include "../../datastructures/datastructures.h"
#include "../Update.h"
class Element;
class Elements;
class Gauss;
class Input;
/*}}}*/

class Material: public Object,public Update{

	public: 
		virtual ~Material(){};
		/*WARNING: input should not be public but it is an easy way to update B from T (using UpdateFromSolution) from Pentas*/

		/*Numerics*/
		virtual void       Configure(Elements* elements)=0;
		virtual Material*  copy2(Element* element)=0;
		virtual IssmDouble GetA()=0;
		virtual IssmDouble GetAbar()=0;
		virtual IssmDouble GetB()=0;
		virtual IssmDouble GetBbar()=0;
		virtual IssmDouble GetD()=0;
		virtual IssmDouble GetDbar()=0;
		virtual IssmDouble GetN()=0;
		virtual void       GetViscosity(IssmDouble* pviscosity,IssmDouble epseff)=0;
		virtual void       GetViscosityBar(IssmDouble* pviscosity,IssmDouble epseff)=0;
		virtual void       GetViscosityComplement(IssmDouble* pviscosity_complement, IssmDouble* pepsilon)=0;
		virtual void       GetViscosityDComplement(IssmDouble* pviscosity_complement, IssmDouble* pepsilon)=0;
		virtual void       GetViscosityDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon)=0;
		virtual void       GetViscosity_B(IssmDouble* pviscosity,IssmDouble epseff)=0;
		virtual void       GetViscosity_D(IssmDouble* pviscosity,IssmDouble epseff)=0;
		virtual void       GetViscosity2dDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon)=0;
		virtual bool       IsDamage()=0;
		virtual bool       IsEnhanced()=0;
		virtual void       ResetHooks()=0;

		virtual void       ViscosityFS(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input)=0;
		virtual void       ViscosityFSDerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon)=0;
		virtual void       ViscosityHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input)=0;
		virtual void       ViscosityHODerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon)=0;
		virtual void       ViscosityL1L2(IssmDouble* pviscosity,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* surf)=0;
		virtual void       ViscositySSA(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input)=0;
		virtual void       ViscositySSADerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon)=0;
		virtual void       ViscosityBFS(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input)=0;
		virtual void       ViscosityBHO(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input)=0;
		virtual void       ViscosityBSSA(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input)=0;

};
#endif
