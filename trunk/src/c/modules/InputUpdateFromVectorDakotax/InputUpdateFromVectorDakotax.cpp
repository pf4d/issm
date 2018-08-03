/*!\file InputUpdateFromVectorDakotax
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromVectorDakotax.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void InputUpdateFromVectorDakotax(FemModel* femmodel,Vector<IssmDouble>* vector, int name, int type){

	IssmDouble* serial_vector=vector->ToMPISerial();
	InputUpdateFromVectorDakotax(femmodel,serial_vector,name, type);

	/*Free ressources:*/
	xDelete<double>(serial_vector);
}

void InputUpdateFromVectorDakotax(FemModel* femmodel,IssmDouble* vector, int name, int type){

	int i;

	/*Update elements, nodes, loads and materials from inputs: */
	for(i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->InputUpdateFromVectorDakota(vector,name,type);
	}
	for(i=0;i<femmodel->loads->Size();i++){
		Load* load=(Load*)femmodel->loads->GetObjectByOffset(i);
		load->InputUpdateFromVectorDakota(vector,name,type);
	}
	for(i=0;i<femmodel->materials->Size();i++){
		Material* material=(Material*)femmodel->materials->GetObjectByOffset(i);
		material->InputUpdateFromVectorDakota(vector,name,type);
	}
}
