#ifndef _CONTAINER_INPUTS_H_
#define _CONTAINER_INPUTS_H_

/*forward declarations */
class Parameters;
class Input;

#include "../../datastructures/datastructures.h"
#include "../../shared/shared.h"

/*! \brief Declaration of Inputs class.
 *
 * Declaration of Inputs class.  Inputs are vector lists (Containers) of Input objects.
 */ 
class Inputs: public DataSet{

	public:

		/*constructors, destructors*/
		Inputs();
		~Inputs();

		/*numerics*/
		int         AddInput(Input* in_input);
		void        AXPY(int inputy_enum, IssmDouble scalar, int inputx_enum);
		void        ChangeEnum(int enumtype,int new_enumtype);
		void        Configure(Parameters* parameters);
		void        ConstrainMin(int constrain_enum, IssmDouble minimum);
		int         DeleteInput(int enum_type);
		void        DuplicateInput(int original_enum,int new_enum);
		Input*      GetInput(int enum_name);
		void        GetInputAverage(IssmDouble* pvalue, int enum_type);
		void        GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void        GetInputValue(bool* pvalue,int enum_type);
		void        GetInputValue(int* pvalue,int enum_type);
		void        GetInputValue(IssmDouble* pvalue,int enum_type);
		IssmDouble  InfinityNorm(int enumtype);
		IssmDouble  Max(int enumtype);
		IssmDouble  MaxAbs(int enumtype);
		IssmDouble  Min(int enumtype);
		IssmDouble  MinAbs(int enumtype);
		Inputs*     SpawnSegInputs(int index1,int index2);
		Inputs*     SpawnSegInputs(int position);
		Inputs*     SpawnTriaInputs(int position);//TO BE REMOVED (replaced by the other one)
		Inputs*     SpawnTriaInputs(int index1,int index2,int index3);

};

#endif //ifndef _INPUTS_H_
