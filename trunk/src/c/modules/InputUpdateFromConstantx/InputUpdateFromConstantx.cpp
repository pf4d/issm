/*!\file InputUpdateFromConstantx
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromConstantx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void InputUpdateFromConstantx(FemModel* femmodel,bool constant, int name){

	int i;
	if(VerboseModule()) _printf0_("   Input updates from constant\n");

	/*Elements and loads drive the update: */
	for(i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->InputUpdateFromConstant(constant,name);
	}

	for(i=0;i<femmodel->loads->Size();i++){
		Load* load=(Load*)femmodel->loads->GetObjectByOffset(i);
		load->InputUpdateFromConstant(constant,name);
	}

	for(i=0;i<femmodel->materials->Size();i++){
		Material* material=(Material*)femmodel->materials->GetObjectByOffset(i);
		material->InputUpdateFromConstant(constant,name);
	}
}
void InputUpdateFromConstantx(FemModel* femmodel,int constant, int name){

	int i;
	if(VerboseModule()) _printf0_("   Input updates from constant\n");

	/*Elements and loads drive the update: */
	for(i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->InputUpdateFromConstant(constant,name);
	}

	for(i=0;i<femmodel->loads->Size();i++){
		Load* load=(Load*)femmodel->loads->GetObjectByOffset(i);
		load->InputUpdateFromConstant(constant,name);
	}

	for(i=0;i<femmodel->materials->Size();i++){
		Material* material=(Material*)femmodel->materials->GetObjectByOffset(i);
		material->InputUpdateFromConstant(constant,name);
	}
}
void InputUpdateFromConstantx(FemModel* femmodel,IssmDouble constant, int name){

	int i;
	if(VerboseModule()) _printf0_("   Input updates from constant\n");

	/*Elements and loads drive the update: */
	for(i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->InputUpdateFromConstant(constant,name);
	}

	for(i=0;i<femmodel->loads->Size();i++){
		Load* load=(Load*)femmodel->loads->GetObjectByOffset(i);
		load->InputUpdateFromConstant(constant,name);
	}

	for(i=0;i<femmodel->materials->Size();i++){
		Material* material=(Material*)femmodel->materials->GetObjectByOffset(i);
		material->InputUpdateFromConstant(constant,name);
	}

}
void InputUpdateFromConstantx(Elements* elements,IssmDouble constant, int name){

	int i;
	if(VerboseModule()) _printf0_("   Input updates from constant\n");

	/*Elements and loads drive the update: */
	for(i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->InputUpdateFromConstant(constant,name);
	}
}
