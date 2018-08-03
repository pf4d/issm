/*!\file ControlInputSetGradientx
 * \brief retrieve gradient from inputs in elements
 */

#include "./ControlInputSetGradientx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void ControlInputSetGradientx(Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters,IssmDouble* gradient){

	/*Intermediaries*/
	int  num_controls;
	int *control_type = NULL;

	/*Retrieve some parameters*/
	parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);

	for(int i=0;i<num_controls;i++){
		for(int j=0;j<elements->Size();j++){
			Element* element=(Element*)elements->GetObjectByOffset(j);
			element->ControlInputSetGradient(gradient,control_type[i],i);
		}
	}

	/*Clean up and return*/
	xDelete<int>(control_type);

}
void ControlInputSetGradientx(Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters,Vector<IssmDouble>* gradient){

	/*Serialize gradient*/
	IssmDouble* serial_gradient=NULL;
	serial_gradient=gradient->ToMPISerial();

	ControlInputSetGradientx(elements,nodes,vertices, loads, materials, parameters,serial_gradient);

	/*Clean up and return*/
	xDelete<IssmDouble>(serial_gradient);
}
