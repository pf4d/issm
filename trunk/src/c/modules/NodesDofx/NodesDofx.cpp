/*!\file NodesDofx
 * \brief: establish degrees of freedom for all nodes
 */

#include "./NodesDofx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void NodesDofx(Nodes* nodes, Parameters* parameters,int configuration_type){

	/*Do we have any nodes for this analysis type? :*/
	if(!nodes->NumberOfNodes(configuration_type)) return;

	/*Do we really need to update dof indexings*/
	if(!nodes->RequiresDofReindexing(configuration_type)) return;

	if(VerboseModule()) _printf0_("   Renumbering degrees of freedom\n");

	/*Ensure that only for each cpu, the partition border nodes only will be taken into account once 
	 * across the cluster. To do so, we flag all the clone nodes: */
	nodes->FlagClones(configuration_type);

	/*Go through all nodes, and build degree of freedom lists. Each node gets a fixed number of dofs. When 
	 *a  node has already been distributed dofs on one cpu, all other cpus with the same node cannot distribute it 
	 *anymore. Use clone field to be sure of that: */
	nodes->DistributeDofs(configuration_type,GsetEnum);
	nodes->DistributeDofs(configuration_type,FsetEnum);
	nodes->DistributeDofs(configuration_type,SsetEnum);

}
