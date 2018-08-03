/*!\file FloatingiceMeltingRatex
 * \brief: calculates Floating ice melting rate
 */

#include "./FloatingiceMeltingRatex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void FloatingiceMeltingRatex(FemModel* femmodel){/*{{{*/

	/*Intermediaties*/
	int  basalforcing_model;

	/*First, get SMB model from parameters*/
	femmodel->parameters->FindParam(&basalforcing_model,BasalforcingsEnum);

	/*branch to correct module*/
	switch(basalforcing_model){
		case FloatingMeltRateEnum:
		case MantlePlumeGeothermalFluxEnum:
			/*Nothing to be done*/
			break;
		case LinearFloatingMeltRateEnum:
			if(VerboseSolution())_printf0_("	call Linear Floating melting rate module\n");
			LinearFloatingiceMeltingRatex(femmodel);
			break;
		case MismipFloatingMeltRateEnum:
			if(VerboseSolution())_printf0_("	call Mismip Floating melting rate module\n");
			MismipFloatingiceMeltingRatex(femmodel);
			break;
		default:
			_error_("Basal forcing model "<<EnumToStringx(basalforcing_model)<<" not supported yet");
	}

}/*}}}*/

void LinearFloatingiceMeltingRatex(FemModel* femmodel){/*{{{*/

	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->LinearFloatingiceMeltingRate();
	}

}/*}}}*/
void MismipFloatingiceMeltingRatex(FemModel* femmodel){/*{{{*/

	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->MismipFloatingiceMeltingRate();
	}

}/*}}}*/
