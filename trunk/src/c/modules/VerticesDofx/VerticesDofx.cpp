/*!\file VerticesDofx
 * \brief: establish degrees of freedom for all vertices: */

#include "./VerticesDofx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void VerticesDofx( Vertices* vertices, Parameters* parameters) {

	/*intermediary: */
	int  numberofvertices;

	/*figure out how many vertices we have: */
	numberofvertices=vertices->NumberOfVertices();

	/*Ensure that only for each cpu, the partition border vertices only will be
	 * taken into account once across the cluster. To do so, we flag all the
	 * clone vertices: */
	vertices->FlagClones(numberofvertices);

	/*Go through all vertices and distribute pids*/
	vertices->DistributePids(numberofvertices); 

}
