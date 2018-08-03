/*!\file Calvingx
 * \brief: compute inverse method gradient
 */

#include "./Calvingx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void Calvingx(FemModel* femmodel){

	/*Recover Calving law Enum*/
	int calvinglaw;
	femmodel->parameters->FindParam(&calvinglaw,CalvingLawEnum);

	/*Calculate calving rate*/
	switch(calvinglaw){
		case DefaultCalvingEnum:
		case CalvingMinthicknessEnum:
			break;
		case CalvingLevermannEnum:
			if(VerboseModule()) _printf0_("   computing Levermann's calving rate\n");
			femmodel->StrainRateparallelx();
			femmodel->StrainRateperpendicularx();
			femmodel->CalvingRateLevermannx();
			break;
		case CalvingDevEnum:
			femmodel->CalvingRateDevx();
			femmodel->ElementOperationx(&Element::CalvingRateDev);
			break;
		default:
			_error_("Caving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
	}
}
