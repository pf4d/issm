/*!\file InputDepthAverageAtBasex
 * \brief: extrude input
 */

#include "./InputDepthAverageAtBasex.h"
#include "../../shared/shared.h"
#include "../../classes/classes.h"
#include "../../toolkits/toolkits.h"

void InputDepthAverageAtBasex(FemModel* femmodel,int original_enum, int new_enum){
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->InputDepthAverageAtBase(original_enum,new_enum);
	}
}
