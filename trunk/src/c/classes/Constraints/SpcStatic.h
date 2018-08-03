/*!\file SpcStatic.h
 * \brief: header file for spc object
 */

#ifndef _SPCStatic_H_
#define _SPCStatic_H_

/*Headers:*/
/*{{{*/
#include "../../datastructures/datastructures.h"
/*}}}*/

class SpcStatic: public Constraint{

	private: 
		int        sid;             /*! id, to track it */
		int        nodeid;          /*!node id          */
		int        dof;             /*!component        */
		IssmDouble value;           /*value             */
		int        analysis_type;
		bool       penalty;         /*Is this a penalty constraint */

	public:

		/*SpcStatic constructors, destructors:{{{*/
		SpcStatic();
		SpcStatic(int sid,int nodeid, int dof,IssmDouble value,int analysis_type);
		~SpcStatic();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int   ObjectEnum();
		/*}}}*/
		/*Constraint virtual functions definitions: {{{*/
		void ActivatePenaltyMethod(void);
		void   ConstrainNode(Nodes* nodes,Parameters* parameters);
		bool   InAnalysis(int analysis_type);
		void   PenaltyDofAndValue(int* dof,IssmDouble* value,Nodes* nodes,Parameters* parameters){_error_("not implemented yet");};
		/*}}}*/
		/*SpcStatic management:{{{ */
		int    GetDof();
		int    GetNodeId();
		IssmDouble GetValue();
		/*}}}*/

};

#endif  /* _SPCStatic_H_*/
