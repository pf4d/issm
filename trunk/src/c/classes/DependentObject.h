/*!\file: DependentObject.h
 * \brief prototype for DependentObject.h
 */ 

#ifndef _DEPENDENTOBJECT_H_
#define  _DEPENDENTOBJECT_H_

/*{{{*/
#include "../datastructures/datastructures.h"
#include "../shared/shared.h"
/*}}}*/

class FemModel;

class DependentObject: public Object{

	public:

		char* name;
		int type;  /*0: scalar, 1: vertex*/
		int index;  /*0: scalar, 1: vertex*/

		/*DependentObject constructors, destructors {{{*/
		DependentObject();
		DependentObject(char* name, int type,int index);
		~DependentObject();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy(void);
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		int   ObjectEnum();
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){/*{{{*/
			_error_("not implemented yet!"); 
		} 
		/*}}}*/
		/*}}}*/

		/*DependentObject methods: */
		int  NumDependents(void);
		void Responsex(IssmDouble* poutput_value,FemModel* femmodel);

};
#endif //ifndef _DEPENDENTOBJECT_H_
