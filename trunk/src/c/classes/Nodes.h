#ifndef _CONTAINER_NODES_H_
#define  _CONTAINER_NODES_H_

#include "../datastructures/datastructures.h"
class Parameters;
class Elements;
class Vertices;
class Loads;
class Nodes;
class Materials;

/*!\brief Declaration of Nodes class.
 *
 * Declaration of Nodes class.  Nodes are vector lists of objects (Containers) of Node objects.
 * Node objects are the degrees of freedom (DOFs) for a particular analysis type (not to be 
 * confused with a vertex, which defines the (x,y,z) location of a point).
 */ 
class Nodes: public DataSet{

	public:

		/*constructors, destructors*/
		Nodes();
		~Nodes();

		/*numerics*/
		void  DistributeDofs(int analysis_type,int SETENUM);
		void  FlagClones(int analysis_type);
		int   MaximumId(void);
		int   MaxNumDofs(int analysis_type,int setenum);
		int   NumberOfDofs(int analysis_type,int setenum);
		int   NumberOfDofsLocal(int analysis_type,int setenum);
		int   NumberOfNodes(int analysis_type);
		int   NumberOfNodes(void);
		void  Ranks(int* ranks,int analysis_type);
		bool  RequiresDofReindexing(int analysis_type);

};

#endif //ifndef _NODES_H_
