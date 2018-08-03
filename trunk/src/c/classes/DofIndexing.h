/*!\file: DofIndexing.h
 * \brief prototype for DofIndexing.h
 */ 

#ifndef _DOFINDEXING_H_
#define _DOFINDEXING_H_

#include "../shared/Numerics/types.h"

class DofIndexing{

	public:

		/*sizes: */
		int gsize;   //number of dofs for a node
		int fsize;   //number of dofs solver for
		int ssize;   //number of constrained dofs

		/*partitioning: */
		bool clone;  //this node is replicated from another one
		bool active; //Is this node active or inactive (all dofs are constrained)
		bool freeze; //this is required for 2d solutions, we never activate nodes that are not on base

		/*boundary conditions sets: */
		bool       *f_set;     //is dof on f-set (on which we solve)
		bool       *s_set;     //is dof on s-set (on which boundary conditions -dirichlet- are applied)
		IssmDouble *svalues;   //list of constraint values. size g_size, for ease of use.

		/*types of dofs: */
		int        *doftype;   //approximation type of the dofs (used only for coupling), size g_size

		/*list of degrees of freedom: */
		int *gdoflist;   //dof list in g_set
		int *fdoflist;   //dof list in f_set
		int *sdoflist;   //dof list in s_set

		/*DofIndexing constructors, destructors {{{*/
		DofIndexing();
		DofIndexing(int g_size);
		void Init(int g_size,int* doftype);
		void InitSet(int setenum);
		DofIndexing(DofIndexing* properties);
		~DofIndexing();
		DofIndexing operator=(const DofIndexing& in);
		/*}}}*/
		/*Object like functionality: {{{*/
		void  copy(const DofIndexing& in);
		void  DeepEcho(void); 
		void  Echo(void); 
		void    Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		/*}}}*/
		/*DofIndexing management: {{{*/
		void Activate(void);
		void Deactivate(void);
		DofIndexing* Spawn(int* indices, int numindices);
		/*}}}*/

};
#endif //ifndef _DOFINDEXING_H_
