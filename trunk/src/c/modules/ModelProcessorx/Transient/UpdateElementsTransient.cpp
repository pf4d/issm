/*
 * UpdateElementsTransient:
 */

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"
#include "../../MeshPartitionx/MeshPartitionx.h"
#include "../ModelProcessorx.h"

void	UpdateElementsTransient(Elements* elements, Parameters* parameters,IoModel* iomodel,int analysis_type){

	bool isgroundingline;
	parameters->FindParam(&isgroundingline,TransientIsgroundinglineEnum);

	if(isgroundingline){
		iomodel->FetchDataToInput(elements,"md.geometry.bed",BedEnum);
	}
}
