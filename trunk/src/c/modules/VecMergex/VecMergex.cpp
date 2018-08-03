/*!\file VecMergex
 * \brief: merge one vector into another
 */

#include "./VecMergex.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
void VecMergex(Vector<IssmDouble>* ug, Vector<IssmDouble>* uf, Nodes* nodes, Parameters* parameters, int SetEnum){

	/*variables: */
	int i;
	int configuration_type;
	IssmDouble* uf_serial=NULL;

	/*retrieve parameters: */
	parameters->FindParam(&configuration_type,ConfigurationTypeEnum);

	/*serialize uf: */
	uf_serial=uf->ToMPISerial();

	/*Do we have any nodes for this configuration? :*/
	if(nodes->NumberOfNodes(configuration_type)){ 

		/*yes. Go through all nodes, and ask them to retrieve values from uf, and plug them into ug: */
		for(i=0;i<nodes->Size();i++){

			Node* node=(Node*)nodes->GetObjectByOffset(i);

			/*Check that this node corresponds to our configuration currently being carried out: */
			if (node->InAnalysis(configuration_type)){

				/*For this object, merge values for enum set SetEnum: */
				node->VecMerge(ug,uf_serial,SetEnum);
			}
		}
	}
	/*Free ressources:*/
	xDelete<IssmDouble>(uf_serial);

	/*Assemble vector: */
	ug->Assemble();
}
