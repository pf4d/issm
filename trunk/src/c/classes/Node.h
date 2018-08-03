/*!\file Node.h
 * \brief: header file for node object
 */

#ifndef _NODE_H_
#define _NODE_H_

/*Headers:*/
/*{{{*/
#include "../datastructures/datastructures.h"
#include "../shared/shared.h"
#include "./DofIndexing.h"
#include "./Update.h"
class  Inputs;
class  Hook;
class  IoModel;
class  DataSet;
class  Vertices;
template <class doubletype> class  Vector;
template <class doubletype> class  Matrix;
class ElementVector;
class ElementMatrix;
/*}}}*/

class Node: public Object{

	private:
		int approximation; //For ice flow models, we need to know what ice flow approximation is employed on this node

	public: 

		int id;    // unique arbitrary id.
		int sid;   // "serial" id (rank of this node if the dataset was serial on 1 cpu)
		int lid;   // "local"  id (rank of this node in current partition)

		int          analysis_enum;
		IssmDouble   coord_system[3][3];
		bool         indexingupdate;
		DofIndexing  indexing;

		/*Node constructors, destructors*/
		Node();
		Node(int node_id,int node_sid,int node_lid,int io_index, IoModel* iomodel,int analysis_enum,int approximation_in);
		~Node();

		/*Object virtual functions definitions:*/
		Object *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int     ObjectEnum();

		/*Node numerical routines*/
		void  Activate(void);
		void  ApplyConstraint(int dof,IssmDouble value);
		void  CreateNodalConstraints(Vector<IssmDouble>* ys);
		void  Deactivate(void);
		void  DistributeDofs(int* pdofcount,int setenum);
		void  DofInFSet(int dof);
		void  DofInSSet(int dof);
		void  FreezeDof(int dof);
		int   GetApproximation();
		void  GetCoordinateSystem(IssmDouble* coord_system_out);
		int   GetDof(int dofindex,int setenum);
		void  GetDofList(int* poutdoflist,int approximation_enum,int setenum);
		void  GetLocalDofList(int* poutdoflist,int approximation_enum,int setenum);
		int   GetNumberOfDofs(int approximation_enum,int setenum);
		void  HardDeactivate(void);
		bool  InAnalysis(int analysis_enum);
		bool  IsActive(void);
		int   IsClone();
		int   IsFloating();
		int   IsGrounded();
		int   Lid(void); 
		void  OffsetDofs(int dofcount,int setenum);
		void  ReindexingDone(void);
		void  RelaxConstraint(int dof);
		bool  RequiresDofReindexing(void);
		void  SetClone(int* minranks);
		void  SetCurrentConfiguration(DataSet* nodes,Vertices* vertices);
		void  ShowTrueDofs(int* truerows,int ncols,int setenum);
		int   Sid(void); 
		void  UpdateCloneDofs(int* alltruerows,int ncols,int setenum);
		void  VecMerge(Vector<IssmDouble>* ug, IssmDouble* vector_serial,int setenum);
		void  VecReduce(Vector<IssmDouble>* vector, IssmDouble* ug_serial,int setnum);
		void  SetApproximation(int in_approximation);
};

/*Methods inherent to Node: */
int* GetGlobalDofList(Node** nodes,int numnodes,int setenum,int approximation);
int* GetLocalDofList(Node** nodes,int numnodes,int setenum,int approximation);
int  GetNumberOfDofs(Node** nodes,int numnodes,int setenum,int approximation);

#endif  /* _NODE_H_ */
