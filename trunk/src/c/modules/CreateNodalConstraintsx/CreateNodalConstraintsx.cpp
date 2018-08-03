/*!\file CreateNodalConstraintsx
 * \brief: establish degrees of freedom for all nodes, and return partitioning vector. Do only once.
 */

#include "./CreateNodalConstraintsx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void CreateNodalConstraintsx( Vector<IssmDouble>** pys, Nodes* nodes,int configuration_type){

	bool  oldalloc  = false;

	/*output: */
	Vector<IssmDouble>* ys=NULL;
	int ssize;
	int slocalsize;

	if(VerboseModule()) _printf0_("   Create nodal constraints\n");

	/*figure out how many dofs we have: */
	ssize=nodes->NumberOfDofs(configuration_type,SsetEnum);
	slocalsize = nodes->NumberOfDofsLocal(configuration_type,SsetEnum);

	/*allocate:*/
	if(oldalloc)
	 ys=new Vector<IssmDouble>(ssize);
	else
	 ys=new Vector<IssmDouble>(slocalsize,ssize);

	/*go through all nodes, and for the ones corresponding to this configuration_type, fill the 
	 * constraints vector with the constraint values: */
	for(int i=0;i<nodes->Size();i++){
		Node* node=(Node*)nodes->GetObjectByOffset(i);
		if (node->InAnalysis(configuration_type)){
			node->CreateNodalConstraints(ys);
		}
	}

	/*Assemble: */
	ys->Assemble();

	/*Assign output pointers: */
	*pys=ys;
}
