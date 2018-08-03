/*!\file Matpar.c
 * \brief: implementation of the Matpar object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*Matpar constructors and destructor*/
Matpar::Matpar(){/*{{{*/
	return;
}
/*}}}*/
Matpar::Matpar(int matpar_mid, IoModel* iomodel){/*{{{*/

	rho_ice                   = 0;
	rho_water                 = 0;
	rho_freshwater            = 0;
	mu_water                  = 0;
	heatcapacity              = 0;
	thermalconductivity       = 0;
	temperateiceconductivity  = 0;
	latentheat                = 0;
	beta                      = 0;
	meltingpoint              = 0;
	referencetemperature      = 0;
	mixed_layer_capacity      = 0;
	thermal_exchange_velocity = 0;
	g                         = 0;
	omega                     = 0;
	desfac                    = 0;
	rlaps                     = 0;
	rlapslgm                  = 0;
	dpermil                   = 0;

	albedo_snow               = 0;
	albedo_ice                = 0;

	sediment_compressibility  = 0;
	sediment_porosity         = 0;
	sediment_thickness        = 0;
	water_compressibility     = 0;

	epl_compressibility       = 0;
	epl_porosity              = 0;
	epl_init_thickness        = 0;
	epl_colapse_thickness     = 0;
	epl_max_thickness         = 0;
	epl_conductivity          = 0;

	lithosphere_shear_modulus = 0;
	lithosphere_density       = 0;
	mantle_shear_modulus      = 0;
	mantle_density            = 0;
	
	earth_density             = 0;

	poisson                   = 0;
	young_modulus             = 0;
	ridging_exponent          = 0;
	cohesion                  = 0;
	internal_friction_coef    = 0;
	compression_coef          = 0;
	traction_coef             = 0;
	time_relaxation_stress    = 0;
	time_relaxation_damage    = 0;

	bool isefficientlayer;
	int  hydrology_model,smb_model,materials_type;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	iomodel->FindConstant(&smb_model,"md.smb.model");
	iomodel->FindConstant(&materials_type,"md.materials.type");

	this->mid = matpar_mid;

	switch(materials_type){
		case MaticeEnum:
		case MatdamageiceEnum:
		case MatenhancediceEnum:
		case MatestarEnum:
			iomodel->FindConstant(&this->rho_ice,"md.materials.rho_ice");
			iomodel->FindConstant(&this->rho_water,"md.materials.rho_water");
			iomodel->FindConstant(&this->rho_freshwater,"md.materials.rho_freshwater");
			iomodel->FindConstant(&this->mu_water,"md.materials.mu_water");
			iomodel->FindConstant(&this->heatcapacity,"md.materials.heatcapacity");
			iomodel->FindConstant(&this->thermalconductivity,"md.materials.thermalconductivity");
			iomodel->FindConstant(&this->temperateiceconductivity,"md.materials.temperateiceconductivity");
			iomodel->FindConstant(&this->latentheat,"md.materials.latentheat");
			iomodel->FindConstant(&this->beta,"md.materials.beta");
			iomodel->FindConstant(&this->meltingpoint,"md.materials.meltingpoint");
			iomodel->FindConstant(&this->referencetemperature,"md.constants.referencetemperature");
			iomodel->FindConstant(&this->mixed_layer_capacity,"md.materials.mixed_layer_capacity");
			iomodel->FindConstant(&this->thermal_exchange_velocity,"md.materials.thermal_exchange_velocity");
			iomodel->FindConstant(&this->g,"md.constants.g");

			switch(smb_model){
				case SMBforcingEnum:
					/*Nothing to add*/
					break;
				case SMBgembEnum:
					iomodel->FindConstant(&this->albedo_ice,"md.smb.aIce");
					iomodel->FindConstant(&this->albedo_snow,"md.smb.aSnow");
					break;
				case SMBpddEnum:
					iomodel->FindConstant(&this->desfac,"md.smb.desfac");
					iomodel->FindConstant(&this->rlaps,"md.smb.rlaps");
					iomodel->FindConstant(&this->rlapslgm,"md.smb.rlapslgm");
					break;
				case SMBd18opddEnum:
					iomodel->FindConstant(&this->desfac,"md.smb.desfac");
					iomodel->FindConstant(&this->rlaps,"md.smb.rlaps");
					iomodel->FindConstant(&this->rlapslgm,"md.smb.rlapslgm");
					iomodel->FindConstant(&this->dpermil,"md.smb.dpermil");					
				case SMBgradientsEnum:
					/*Nothing to add*/
					break;
				case SMBgradientselaEnum:
					/*Nothing to add*/
					break;
				case SMBhenningEnum:
					/*Nothing to add*/
					break;
				case SMBcomponentsEnum:
					/*Nothing to add*/
					break;
				case SMBmeltcomponentsEnum:
					/*Nothing to add*/
					break;
				default:
					_error_("Surface mass balance model "<<EnumToStringx(smb_model)<<" not supported yet");
			}
			if(hydrology_model==HydrologydcEnum){
				iomodel->FindConstant(&this->sediment_compressibility,"md.hydrology.sediment_compressibility");
				iomodel->FindConstant(&this->sediment_porosity,"md.hydrology.sediment_porosity");
				iomodel->FindConstant(&this->sediment_thickness,"md.hydrology.sediment_thickness");
				iomodel->FindConstant(&this->water_compressibility,"md.hydrology.water_compressibility");
				iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");

				if(isefficientlayer){
					iomodel->FindConstant(&this->epl_compressibility,"md.hydrology.epl_compressibility");
					iomodel->FindConstant(&this->epl_porosity,"md.hydrology.epl_porosity");
					iomodel->FindConstant(&this->epl_init_thickness,"md.hydrology.epl_initial_thickness");
					iomodel->FindConstant(&this->epl_colapse_thickness,"md.hydrology.epl_colapse_thickness");
					iomodel->FindConstant(&this->epl_max_thickness,"md.hydrology.epl_max_thickness");
					iomodel->FindConstant(&this->epl_conductivity,"md.hydrology.epl_conductivity");
				}
			}
			else if(hydrology_model==HydrologyshreveEnum){
				/*Nothing to add*/
			}
			else if(hydrology_model==HydrologysommersEnum){
				/*Nothing to add*/
			}
			else{
				_error_("Hydrology model "<<EnumToStringx(hydrology_model)<<" not supported yet");
			}

			/*gia: */
			iomodel->FindConstant(&this->lithosphere_shear_modulus,"md.materials.lithosphere_shear_modulus");
			iomodel->FindConstant(&this->lithosphere_density,"md.materials.lithosphere_density");
			iomodel->FindConstant(&this->mantle_shear_modulus,"md.materials.mantle_shear_modulus");
			iomodel->FindConstant(&this->mantle_density,"md.materials.mantle_density");

			/*slr:*/
			iomodel->FindConstant(&this->earth_density,"md.materials.earth_density");

			break;
		default:
			_error_("Material "<< EnumToStringx(materials_type) <<" not supported yet");
	}
}
/*}}}*/
Matpar::~Matpar(){/*{{{*/
	return;
}
/*}}}*/
void Matpar::SetMid(int matpar_mid){/*{{{*/
	this->mid=matpar_mid;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Matpar::copy() {/*{{{*/

	/*Output*/
	Matpar* matpar;

	/*Initialize output*/
	matpar=new Matpar(*this);

	/*copy fields: */
	matpar->mid=this->mid;
	matpar->rho_ice=this->rho_ice;
	matpar->rho_water=this->rho_water;
	matpar->rho_freshwater=this->rho_freshwater;
	matpar->mu_water=this->mu_water;
	matpar->heatcapacity=this->heatcapacity;
	matpar->thermalconductivity=this->thermalconductivity;
	matpar->temperateiceconductivity=this->temperateiceconductivity;
	matpar->latentheat=this->latentheat;
	matpar->beta=this->beta;
	matpar->meltingpoint=this->meltingpoint;
	matpar->referencetemperature=this->referencetemperature;
	matpar->mixed_layer_capacity=this->mixed_layer_capacity;
	matpar->thermal_exchange_velocity=this->thermal_exchange_velocity;
	matpar->g=this->g;
	matpar->desfac=this->desfac;
	matpar->rlaps=this->rlaps;
	matpar->rlapslgm=this->rlapslgm;
	matpar->dpermil=this->dpermil;

	matpar->sediment_compressibility=this->sediment_compressibility;
	matpar->sediment_porosity=this->sediment_porosity;
	matpar->sediment_thickness=this->sediment_thickness;
	matpar->water_compressibility=this->water_compressibility;

	matpar->epl_compressibility=this->epl_compressibility;
	matpar->epl_porosity=this->epl_porosity;
	matpar->epl_init_thickness=this->epl_init_thickness;
	matpar->epl_colapse_thickness=this->epl_colapse_thickness;
	matpar->epl_max_thickness=this->epl_max_thickness;
	matpar->epl_conductivity=this->epl_conductivity;

	matpar->lithosphere_shear_modulus=this->lithosphere_shear_modulus;
	matpar->lithosphere_density=this->lithosphere_density;
	matpar->mantle_shear_modulus=this->mantle_shear_modulus;
	matpar->mantle_density=this->mantle_density;
	
	matpar->earth_density=this->earth_density;

	return matpar;
}
/*}}}*/
void Matpar::DeepEcho(void){/*{{{*/

	this->Echo();
}		
/*}}}*/
void Matpar::Echo(void){/*{{{*/

	_printf_("Matpar:\n");
	_printf_("   mid: " << mid << "\n");
	_printf_("   rho_ice: " << rho_ice << "\n");
	_printf_("   rho_water: " << rho_water << "\n");
	_printf_("   rho_freshwater: " << rho_freshwater << "\n");
	_printf_("   mu_water: " << mu_water << "\n");
	_printf_("   heatcapacity: " << heatcapacity << "\n");
	_printf_("   thermalconductivity: " << thermalconductivity << "\n");
	_printf_("   temperateiceconductivity: " << temperateiceconductivity << "\n");
	_printf_("   latentheat: " << latentheat << "\n");
	_printf_("   beta: " << beta << "\n");
	_printf_("   meltingpoint: " << meltingpoint << "\n");
	_printf_("   referencetemperature: " << referencetemperature << "\n");
	_printf_("   mixed_layer_capacity: " << mixed_layer_capacity << "\n");
	_printf_("   thermal_exchange_velocity: " << thermal_exchange_velocity << "\n");
	_printf_("   g: " << g << "\n");
	_printf_("   omega: " << omega << "\n");
	_printf_("   desfac: " << desfac << "\n");
	_printf_("   rlaps: " << rlaps << "\n");
	_printf_("   rlapslgm: " << rlapslgm << "\n");
	_printf_("   dpermil: " << dpermil << "\n");
	_printf_("   albedo_ice: " << albedo_ice << "\n");
	_printf_("   albedo_snow: " << albedo_snow << "\n");
	_printf_("   sediment_compressibility: " << sediment_compressibility << "\n");
	_printf_("   sediment_porosity: " << sediment_porosity << "\n");
	_printf_("   sediment_thickness: " << sediment_thickness << "\n");
	_printf_("   water_compressibility: " << water_compressibility << "\n");	
	_printf_("   epl_compressibility: " << epl_compressibility << "\n");	
	_printf_("   epl_porosity: " << epl_porosity << "\n");	
	_printf_("   epl_init_thickness: " << epl_init_thickness << "\n");	
	_printf_("   epl_colapse_thickness: " << epl_colapse_thickness << "\n");	
	_printf_("   epl_max_thickness: " << epl_max_thickness << "\n");	
	_printf_("   epl_conductivity: " << epl_conductivity << "\n");	
	_printf_("   lithosphere_shear_modulus: " << lithosphere_shear_modulus << "\n");	
	_printf_("   lithosphere_density: " << lithosphere_density << "\n");	
	_printf_("   mantle_shear_modulus: " << mantle_shear_modulus << "\n");	
	_printf_("   mantle_density: " << mantle_density << "\n");	
	_printf_("   earth_density: " << earth_density << "\n");	
	_printf_("   poisson: " << poisson << "\n");	
	_printf_("   young_modulus: " << young_modulus << "\n");	
	_printf_("   ridging_exponent: " << ridging_exponent << "\n");	
	_printf_("   cohesion: " << cohesion << "\n");	
	_printf_("   internal_friction_coef: " << internal_friction_coef << "\n");	
	_printf_("   compression_coef: " << compression_coef << "\n");	
	_printf_("   traction_coef: " << traction_coef << "\n");	
	_printf_("   time_relaxation_stress: " << time_relaxation_stress << "\n");	
	_printf_("   time_relaxation_damage: " << time_relaxation_damage << "\n");	
	return;
}
/*}}}*/
int  Matpar::Id(void){ return mid; }/*{{{*/
/*}}}*/
void Matpar::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(MatparEnum);

	MARSHALLING(mid);
	MARSHALLING(rho_ice);
	MARSHALLING(rho_water);
	MARSHALLING(rho_freshwater);
	MARSHALLING(mu_water);
	MARSHALLING(heatcapacity);
	MARSHALLING(thermalconductivity);
	MARSHALLING(temperateiceconductivity);
	MARSHALLING(latentheat);
	MARSHALLING(beta);
	MARSHALLING(meltingpoint);
	MARSHALLING(referencetemperature);
	MARSHALLING(mixed_layer_capacity);
	MARSHALLING(thermal_exchange_velocity);
	MARSHALLING(g);
	MARSHALLING(omega);
	MARSHALLING(desfac);
	MARSHALLING(rlaps);
	MARSHALLING(rlapslgm);
	MARSHALLING(dpermil);

	//hydrology Dual Porous Continuum:
	MARSHALLING(sediment_compressibility);
	MARSHALLING(sediment_porosity);
	MARSHALLING(sediment_thickness);
	MARSHALLING(water_compressibility);

	MARSHALLING(epl_compressibility);
	MARSHALLING(epl_porosity);
	MARSHALLING(epl_init_thickness);
	MARSHALLING(epl_colapse_thickness);
	MARSHALLING(epl_max_thickness);
	MARSHALLING(epl_conductivity);

	//gia:
	MARSHALLING(lithosphere_shear_modulus);
	MARSHALLING(lithosphere_density);
	MARSHALLING(mantle_shear_modulus);
	MARSHALLING(mantle_density);
	
	//slr:
	MARSHALLING(earth_density);

	//Sea ice:
	MARSHALLING(poisson);
	MARSHALLING(young_modulus);
	MARSHALLING(ridging_exponent);
	MARSHALLING(cohesion);
	MARSHALLING(internal_friction_coef);
	MARSHALLING(compression_coef);
	MARSHALLING(traction_coef);
	MARSHALLING(time_relaxation_stress);
	MARSHALLING(time_relaxation_damage);

}
/*}}}*/
int  Matpar::ObjectEnum(void){/*{{{*/

	return MatparEnum;

}
/*}}}*/

/*Update virtual functions definitions:*/
void   Matpar::InputUpdateFromConstant(IssmDouble constant, int name){/*{{{*/

	switch(name){
		case MaterialsRhoIceEnum:
			this->rho_ice=constant;
			break;
		case MaterialsRhoSeawaterEnum:
			this->rho_water=constant;
			break;
		case MaterialsRhoFreshwaterEnum:
			this->rho_freshwater=constant;
			break;
		case MaterialsMuWaterEnum:
			this->mu_water=constant;
			break;
		case MaterialsHeatcapacityEnum:
			this->heatcapacity=constant;
			break;
	  	case MaterialsThermalconductivityEnum:
			this->thermalconductivity=constant;
			break;
	  	case MaterialsTemperateiceconductivityEnum:
			this->temperateiceconductivity=constant;
			break;
		case  MaterialsLatentheatEnum:
			this->latentheat=constant;
			break;
		case  MaterialsBetaEnum:
			this->beta=constant;
			break;
		case  MaterialsMeltingpointEnum:
			this->meltingpoint=constant;
			break;
		case  ConstantsReferencetemperatureEnum:
			this->referencetemperature=constant;
			break;
		case  MaterialsMixedLayerCapacityEnum:
			this->mixed_layer_capacity=constant;
			break;
		case  MaterialsThermalExchangeVelocityEnum:
			this->thermalconductivity=constant;
			break;
		case  ConstantsGEnum:
			this->g=constant;
			break;
		case  SmbDesfacEnum:
			this->desfac=constant;
			break;
		case SmbRlapsEnum:
			this->rlaps=constant;
			break;
		case SmbRlapslgmEnum:
			this->rlapslgm=constant;
			break;
		case  SmbDpermilEnum:
			this->dpermil=constant;
			break;
		default: 
			break;
	}

}
/*}}}*/
void   Matpar::InputUpdateFromConstant(int constant, int name){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void   Matpar::InputUpdateFromConstant(bool constant, int name){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Matpar::InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrows, int ncols,int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void   Matpar::InputUpdateFromVector(IssmDouble* vector, int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void   Matpar::InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/

/*Matpar management: */
void       Matpar::Configure(Elements* elementsin){/*{{{*/

	/*nothing done yet!*/

}
/*}}}*/
void       Matpar::EnthalpyToThermal(IssmDouble* ptemperature,IssmDouble* pwaterfraction,IssmDouble enthalpy,IssmDouble pressure){/*{{{*/

	/*Ouput*/
	IssmDouble temperature,waterfraction;

	if(enthalpy<PureIceEnthalpy(pressure)){
		temperature=referencetemperature+enthalpy/heatcapacity;
		waterfraction=0.;
	}
	else{
		temperature=TMeltingPoint(pressure);
		waterfraction=(enthalpy-PureIceEnthalpy(pressure))/latentheat;
	}

	/*Assign output pointers:*/
	*pwaterfraction=waterfraction;
	*ptemperature=temperature;
}
/*}}}*/
IssmDouble Matpar::GetEnthalpyDiffusionParameter(IssmDouble enthalpy,IssmDouble pressure){/*{{{*/
	if (enthalpy<PureIceEnthalpy(pressure))
		return thermalconductivity/heatcapacity;
	else
		return temperateiceconductivity/heatcapacity;
}
/*}}}*/
IssmDouble Matpar::GetEnthalpyDiffusionParameterVolume(int numvertices,IssmDouble* enthalpy,IssmDouble* pressure){/*{{{*/

	int         iv;
	IssmDouble  lambda;                 // fraction of cold ice
	IssmDouble  kappa,kappa_c,kappa_t;  //enthalpy conductivities
	IssmDouble  Hc,Ht;
	IssmDouble* PIE   = xNew<IssmDouble>(numvertices);
	IssmDouble* dHpmp = xNew<IssmDouble>(numvertices);

	for(iv=0; iv<numvertices; iv++){
		PIE[iv]=PureIceEnthalpy(pressure[iv]);
		dHpmp[iv]=enthalpy[iv]-PIE[iv];
	}

	bool allequalsign=true;
	if(dHpmp[0]<0)
		for(iv=1; iv<numvertices;iv++) allequalsign=(allequalsign && (dHpmp[iv]<0));
	else
		for(iv=1; iv<numvertices;iv++) allequalsign=(allequalsign && (dHpmp[iv]>=0));

	if(allequalsign){
		kappa=GetEnthalpyDiffusionParameter(enthalpy[0], pressure[0]);
	}
	else {
		/* return harmonic mean of thermal conductivities, weighted by fraction of cold/temperate ice,
		 cf Patankar 1980, pp44 */
		kappa_c=GetEnthalpyDiffusionParameter(PureIceEnthalpy(0.)-1.,0.);
		kappa_t=GetEnthalpyDiffusionParameter(PureIceEnthalpy(0.)+1.,0.);
		Hc=0.; Ht=0.;
		for(iv=0; iv<numvertices;iv++){
			if(enthalpy[iv]<PIE[iv])
			 Hc+=(PIE[iv]-enthalpy[iv]);
			else
			 Ht+=(enthalpy[iv]-PIE[iv]);
		}
		_assert_((Hc+Ht)>0.);
		lambda = Hc/(Hc+Ht);
		kappa  = 1./(lambda/kappa_c + (1.-lambda)/kappa_t);
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(PIE);
	xDelete<IssmDouble>(dHpmp);
	return kappa;
}
/*}}}*/
IssmDouble Matpar::GetMaterialParameter(int enum_in){/*{{{*/

	switch(enum_in){
		case MaterialsRhoIceEnum:                    return this->rho_ice;
		case MaterialsRhoSeawaterEnum:               return this->rho_water;
		case MaterialsRhoFreshwaterEnum:             return this->rho_freshwater;
		case MaterialsMuWaterEnum:                   return this->mu_water;
		case MaterialsHeatcapacityEnum:              return this->heatcapacity;
		case MaterialsThermalconductivityEnum:       return this->thermalconductivity;
		case MaterialsTemperateiceconductivityEnum:  return this->temperateiceconductivity;
		case MaterialsLatentheatEnum:                return this->latentheat;
		case MaterialsBetaEnum:                      return this->beta;
		case MaterialsMeltingpointEnum:              return this->meltingpoint;
		case ConstantsReferencetemperatureEnum:      return this->referencetemperature;
		case MaterialsMixedLayerCapacityEnum:        return this->mixed_layer_capacity;
		case MaterialsThermalExchangeVelocityEnum:   return this->thermal_exchange_velocity;
		case HydrologydcSedimentPorosityEnum:        return this->sediment_porosity;
		case HydrologydcSedimentThicknessEnum:       return this->sediment_thickness;
		case HydrologydcSedimentCompressibilityEnum: return this->sediment_compressibility;
		case HydrologydcEplPorosityEnum:             return this->epl_porosity;
		case HydrologydcEplCompressibilityEnum:      return this->epl_compressibility;
		case HydrologydcEplConductivityEnum:         return this->epl_conductivity;
		case HydrologydcEplInitialThicknessEnum:     return this->epl_init_thickness;
		case HydrologydcEplColapseThicknessEnum:     return this->epl_colapse_thickness;
		case HydrologydcEplMaxThicknessEnum:         return this->epl_max_thickness;
		case HydrologydcWaterCompressibilityEnum:    return this->water_compressibility;
		case ConstantsGEnum:                         return this->g;
		case SmbDesfacEnum:              return this->desfac;
		case SmbRlapsEnum:               return this->rlaps;
		case SmbRlapslgmEnum:            return this->rlapslgm;
		case SmbDpermilEnum:             return this->dpermil;
		case MaterialsLithosphereShearModulusEnum:   return this->lithosphere_shear_modulus;
		case MaterialsLithosphereDensityEnum:        return this->lithosphere_density;
		case MaterialsMantleDensityEnum:             return this->mantle_density;
		case MaterialsMantleShearModulusEnum:        return this->mantle_shear_modulus;
		case MaterialsEarthDensityEnum:              return this->earth_density;
		default: _error_("Enum "<<EnumToStringx(enum_in)<<" not supported yet");
	}

}
/*}}}*/
IssmDouble Matpar::PureIceEnthalpy(IssmDouble pressure){/*{{{*/
	return heatcapacity*(TMeltingPoint(pressure)-referencetemperature);
}
/*}}}*/
void       Matpar::ResetHooks(){/*{{{*/

	//Nothing to be done
	return;
}
/*}}}*/
void       Matpar::ThermalToEnthalpy(IssmDouble * penthalpy,IssmDouble temperature,IssmDouble waterfraction,IssmDouble pressure){/*{{{*/

	/*Ouput*/
	IssmDouble enthalpy;

	if(temperature<TMeltingPoint(pressure)){
		enthalpy=heatcapacity*(temperature-referencetemperature);
	}
	else{
		enthalpy=PureIceEnthalpy(pressure)+latentheat*waterfraction;
	}

	/*Assign output pointers:*/
	*penthalpy=enthalpy;
}
/*}}}*/
IssmDouble Matpar::TMeltingPoint(IssmDouble pressure){/*{{{*/
	return meltingpoint-beta*pressure;
}
/*}}}*/
