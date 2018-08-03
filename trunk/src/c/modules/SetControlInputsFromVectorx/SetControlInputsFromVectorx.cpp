/*!\file SetControlInputsFromVectorx
 * \brief retrieve vector from inputs in elements
 */

#include "./SetControlInputsFromVectorx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void SetControlInputsFromVectorx(FemModel* femmodel,IssmDouble* vector){

	int  num_controls;
	int *control_type = NULL;

	/*Retrieve some parameters*/
	femmodel->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	femmodel->parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);

	for(int i=0;i<num_controls;i++){
		for(int j=0;j<femmodel->elements->Size();j++){
			Element* element=(Element*)femmodel->elements->GetObjectByOffset(j);
			element->SetControlInputsFromVector(vector,control_type[i],i);
		}
	}

	xDelete<int>(control_type);
}

void SetControlInputsFromVectorx(FemModel* femmodel,Vector<IssmDouble>* vector){

	IssmDouble* serial_vector=vector->ToMPISerial();
	SetControlInputsFromVectorx(femmodel,serial_vector);
	xDelete<IssmDouble>(serial_vector);
}
