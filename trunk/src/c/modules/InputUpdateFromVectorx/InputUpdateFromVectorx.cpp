/*!\file InputUpdateFromVectorx
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromVectorx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void InputUpdateFromVectorx(FemModel* femmodel,Vector<IssmDouble>* vector, int name, int type){

	IssmDouble* serial_vector=vector->ToMPISerial();
	InputUpdateFromVectorx(femmodel,serial_vector,name,type);
	xDelete<IssmDouble>(serial_vector);
}

void InputUpdateFromVectorx(FemModel* femmodel,IssmDouble* vector, int name, int type){

	int i;

	/*Update elements, nodes, loads and materials from inputs: */
	for(i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->InputUpdateFromVector(vector,name,type);
	}
	for(i=0;i<femmodel->loads->Size();i++){
		Load* load=(Load*)femmodel->loads->GetObjectByOffset(i);
		load->InputUpdateFromVector(vector,name,type);
	}
	for(i=0;i<femmodel->materials->Size();i++){
		Material* material=(Material*)femmodel->materials->GetObjectByOffset(i);
		material->InputUpdateFromVector(vector,name,type);
	}
}
