/*!\file GeothermalFluxx
 * \brief: calculates Geothermal heat flux 
 */

#include "./GeothermalFluxx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void GeothermalFluxx(FemModel* femmodel){/*{{{*/

	/*Intermediaties*/
	int  basalforcing_model;

	/*First, get SMB model from parameters*/
	femmodel->parameters->FindParam(&basalforcing_model,BasalforcingsEnum);

	/*branch to correct module*/
	switch(basalforcing_model){
		case FloatingMeltRateEnum:
		case MismipFloatingMeltRateEnum:
		case LinearFloatingMeltRateEnum:
			/*Nothing to be done*/
			break;
		case MantlePlumeGeothermalFluxEnum:
			if(VerboseSolution())_printf0_("	call Mantle Plume Geothermal Flux module\n");
			MantlePlumeGeothermalFluxx(femmodel);
			break;
		default:
			_error_("Basal forcing model "<<EnumToStringx(basalforcing_model)<<" not supported yet");
	}

}/*}}}*/

void MantlePlumeGeothermalFluxx(FemModel* femmodel){/*{{{*/

	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->MantlePlumeGeothermalFlux();
	}

}/*}}}*/
