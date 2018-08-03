/*!\file InputExtrudex
 * \brief: extrude input
 */

#include "./InputExtrudex.h"
#include "../../shared/shared.h"
#include "../../classes/classes.h"
#include "../../toolkits/toolkits.h"

void InputExtrudex(FemModel* femmodel,int input_enum,int start){
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->InputExtrude(input_enum,start);
	}
}
