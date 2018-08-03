#ifndef _CONTAINER_LOADS_H_
#define  _CONTAINER_LOADS_H_

/*forward declarations */
#include "../../datastructures/datastructures.h"
class Materials;
class Parameters;
class Elements;
class Vertices;
class Nodes;

/*!\brief Declaration of Loads class.
 *
 * Declaration of Loads class.  Loads are vector lists (Containers) of Load objects.
 */ 
class Loads: public DataSet{

	public:

		/*constructors, destructors*/
		Loads();
		~Loads();

		/*numerics*/
		void  Configure(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters);
		bool  IsPenalty(int analysis);
		int   MaxNumNodes(int analysis);
		int   NumberOfLoads(void);
		int   NumberOfLoads(int analysis);
		void  ResetHooks();
		void  SetCurrentConfiguration(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters);
		int   Size(int analysis);
		int   Size(void);

};

#endif //ifndef _LOADS_H_
