/*!\file Matpar.h
 * \brief: header file for matpar object
 */

#ifndef _MATPAR_H_
#define _MATPAR_H_

/*Headers:*/
/*{{{*/
#include "./Material.h"
class IoModel;
/*}}}*/

class Matpar: public Material{

	private: 
		int	      mid;
		IssmDouble  rho_ice; 
		IssmDouble  rho_water;
		IssmDouble  rho_freshwater;
		IssmDouble  mu_water;
		IssmDouble  heatcapacity;
		IssmDouble  thermalconductivity;
		IssmDouble  temperateiceconductivity;
		IssmDouble  latentheat;
		IssmDouble  beta;
		IssmDouble  meltingpoint;
		IssmDouble  referencetemperature;
		IssmDouble  mixed_layer_capacity;
		IssmDouble  thermal_exchange_velocity;
		IssmDouble  g;
		IssmDouble  omega;
		IssmDouble  desfac;
		IssmDouble  rlaps;
		IssmDouble  rlapslgm;
		IssmDouble  dpermil;

		/*albedo: */
		IssmDouble albedo_ice;
		IssmDouble albedo_snow;

		/*hydrology Dual Porous Continuum: */	 
		IssmDouble  sediment_compressibility;
		IssmDouble  sediment_porosity;	 
		IssmDouble  sediment_thickness;
		IssmDouble  water_compressibility;

		IssmDouble  epl_compressibility;
		IssmDouble  epl_porosity;
		IssmDouble  epl_init_thickness;
		IssmDouble  epl_colapse_thickness;
		IssmDouble  epl_max_thickness;
		IssmDouble  epl_conductivity;	 

		/*gia: */
		IssmDouble lithosphere_shear_modulus;
		IssmDouble lithosphere_density;
		IssmDouble mantle_shear_modulus;
		IssmDouble mantle_density;

		/*slr:*/
		IssmDouble earth_density;

		/*Sea ice*/
		IssmDouble poisson;
		IssmDouble young_modulus;
		IssmDouble ridging_exponent;
		IssmDouble cohesion;
		IssmDouble internal_friction_coef;
		IssmDouble compression_coef;
		IssmDouble traction_coef;
		IssmDouble time_relaxation_stress;
		IssmDouble time_relaxation_damage;

	public:
		Matpar();
		Matpar(int matpar_id, IoModel* iomodel);
		~Matpar();
		void SetMid(int matpar_mid);

		/*Object virtual functions definitions:{{{ */
		Object *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int     ObjectEnum();
		/*}}}*/
		/*Update virtual functions resolution: {{{*/
		void   InputUpdateFromConstant(IssmDouble constant, int name);
		void   InputUpdateFromConstant(int constant, int name);
		void   InputUpdateFromConstant(bool constant, int name);
		void   InputUpdateFromIoModel(int index, IoModel* iomodel){_error_("not implemented");};
		void   InputUpdateFromMatrixDakota(IssmDouble* matrix,int nrows,int ncols, int name, int type);
		void   InputUpdateFromVector(IssmDouble* vector, int name, int type);
		void   InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type);
		/*}}}*/
		/*Material virtual functions resolution: {{{*/
		Material*  copy2(Element* element){_error_("not implemented");};
		void       Configure(Elements* elements);
		void       GetViscosity(IssmDouble* pviscosity,IssmDouble eps_eff){_error_("not supported");};
		void       GetViscosityBar(IssmDouble* pviscosity,IssmDouble eps_eff){_error_("not supported");};
		void       GetViscosityComplement(IssmDouble* pviscosity_complement, IssmDouble* pepsilon){_error_("not supported");};
		void       GetViscosityDComplement(IssmDouble* pviscosity_complement, IssmDouble* pepsilon){_error_("not supported");};
		void       GetViscosityDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon){_error_("not supported");};
		void       GetViscosity_B(IssmDouble* pviscosity,IssmDouble eps_eff){_error_("not supported");};
		void       GetViscosity_D(IssmDouble* pviscosity,IssmDouble eps_eff){_error_("not supported");};
		void       GetViscosity2dDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon){_error_("not supported");};
		IssmDouble GetA(){_error_("not supported");};
		IssmDouble GetAbar(){_error_("not supported");};
		IssmDouble GetB(){_error_("not supported");};
		IssmDouble GetBbar(){_error_("not supported");};
		IssmDouble GetD(){_error_("not supported");};
		IssmDouble GetDbar(){_error_("not supported");};
		IssmDouble GetN(){_error_("not supported");};
		bool       IsDamage(){_error_("not supported");};
		bool       IsEnhanced(){_error_("not supported");};
		void       ResetHooks();

		void       ViscosityFS(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){_error_("not supported");};
		void       ViscosityFSDerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon){_error_("not supported");};
		void       ViscosityHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){_error_("not supported");};
		void       ViscosityHODerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon){_error_("not supported");};
		void       ViscosityL1L2(IssmDouble* pviscosity,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* surf){_error_("not supported");};
		void       ViscositySSA(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){_error_("not supported");};
		void       ViscositySSADerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon){_error_("not supported");};
		void       ViscosityBFS(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){_error_("not supported");};
		void       ViscosityBHO(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){_error_("not supported");};
		void       ViscosityBSSA(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){_error_("not supported");};
		/*}}}*/
		/*Numerics: {{{*/
		void       EnthalpyToThermal(IssmDouble* ptemperature,IssmDouble* pwaterfraction,IssmDouble enthalpy,IssmDouble pressure);
		IssmDouble GetEnthalpyDiffusionParameter(IssmDouble enthalpy,IssmDouble pressure);
		IssmDouble GetEnthalpyDiffusionParameterVolume(int numvertices,IssmDouble* enthalpy,IssmDouble* pressure);
		IssmDouble GetMaterialParameter(int in_enum); 
		IssmDouble PureIceEnthalpy(IssmDouble pressure);
		void       ThermalToEnthalpy(IssmDouble* penthalpy,IssmDouble temperature,IssmDouble waterfraction,IssmDouble pressure);
		IssmDouble TMeltingPoint(IssmDouble pressure);
		/*}}}*/

};

#endif  /* _MATPAR_H_ */
