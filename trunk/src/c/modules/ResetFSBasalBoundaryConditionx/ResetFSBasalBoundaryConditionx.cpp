/*!\file ResetFSBasalBoundaryConditionx
 * \brief: reset coordinate system for full-FS: tangential to the bedrock
 */

#include "./ResetFSBasalBoundaryConditionx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void ResetFSBasalBoundaryConditionx(Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads,Materials* materials,Parameters* parameters){

	Element *element = NULL;

	for (int i=0;i<elements->Size();i++){
		element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->ResetFSBasalBoundaryCondition();
	}

}
