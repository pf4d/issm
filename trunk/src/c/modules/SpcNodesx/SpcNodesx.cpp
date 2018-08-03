/*!\file SpcNodesx
 * \brief: establish single point constraints on all nodes, as well as constraints vector.
 */

#include "./SpcNodesx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void SpcNodesx(Nodes* nodes,Constraints* constraints,Parameters* parameters, int analysis_type){

	for(int i=0;i<constraints->Size();i++){

		Constraint* constraint=(Constraint*)constraints->GetObjectByOffset(i);

		/*Check this constraint belongs to this analysis: */
		if(constraint->InAnalysis(analysis_type)){

			/*Ok, apply constraint onto corresponding node: */
			constraint->ConstrainNode(nodes,parameters);
		}
	}

}
