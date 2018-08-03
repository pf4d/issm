/*!\file InputDuplicatex
 * \brief: duplicte  an input inside the elements, onto another, and wipe it off.
 */

#include "./InputDuplicatex.h"
#include "../../shared/shared.h"
#include "../../classes/classes.h"
#include "../../toolkits/toolkits.h"

void InputDuplicatex(FemModel* femmodel,int original_enum, int new_enum){
	/*Go through elemnets, and ask to reinitialie the input: */
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->InputDuplicate(original_enum,new_enum);
	}
}
