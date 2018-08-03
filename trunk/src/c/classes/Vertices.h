#ifndef _CONTAINER_VERTICES_H_
#define  _CONTAINER_VERTICES_H_

/*forward declarations */
#include "../datastructures/datastructures.h"
#include "../shared/shared.h"

/*!\brief Declaration of Vertices class.
 *
 * Declaration of Vertices class.  Vertices are vector lists (Containers) of Vertex objects.
 * A vertex is a set of (x,y,z) coordinates defining the location of points in the mesh (not
 * to be confused with a node, which is a degree of freedom (DOF) for a particular analysis).
 */ 
class Vertices: public DataSet{

	public:

		/*constructors, destructors:*/ 
		Vertices();
		~Vertices();

		/*numerics:*/
		void  DistributePids(int numberofnodes);
		void  FlagClones(int numberofnodes);
		int   NumberOfVertices(void);
		void  Ranks(int* ranks);
		IssmDouble* ToXYZ(void);
};

#endif //ifndef _VERTICES_H_
