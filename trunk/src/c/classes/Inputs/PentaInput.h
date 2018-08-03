/*! \file PentaInput.h 
 *  \brief: header file for PentaInput object
 */

#ifndef _PENTAINPUT_H_
#define _PENTAINPUT_H_

/*Headers:*/
/*{{{*/
#include "./Input.h"
#include "../Elements/PentaRef.h"
class Gauss;
class Gauss;
/*}}}*/

class PentaInput: public Input, public PentaRef{

	public:
		int         enum_type;
		int         interpolation_type;
		IssmDouble* values;

		/*PentaInput constructors, destructors*/
		PentaInput();
		PentaInput(int enum_type,IssmDouble* values,int element_type_in);
		~PentaInput();

		/*Object virtual functions definitions */
		Object *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int     ObjectEnum();

		/*PentaInput management*/
		int   InstanceEnum();
		Input* PointwiseDivide(Input* inputB);
		Input* PointwiseMin(Input* inputB);
		Input* PointwiseMax(Input* inputB);
		int  GetResultInterpolation(void);
		int  GetResultNumberOfNodes(void);
		int  GetResultArraySize(void){return 1;};
		void ResultToPatch(IssmDouble* values,int nodesperelement,int sid);
		void AddTimeValues(IssmDouble* values,int step,IssmDouble time){_error_("not supported yet");};
		void Configure(Parameters* parameters);
		/*numerics*/
		void AXPY(Input* xinput,IssmDouble scalar);
		void ChangeEnum(int newenumtype);
		void Constrain(IssmDouble cm_min, IssmDouble cm_max);
		void ConstrainMin(IssmDouble minimum);
		void Extrude(int start);
		void GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes){_error_("not implemented yet");};
		void GetInputAverage(IssmDouble* pvalue);
		void GetInputDerivativeAverageValue(IssmDouble* derivativevalues, IssmDouble* xyz_list){_error_("not implemented yet");};
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list,Gauss* gauss);
		void GetInputValue(bool* pvalue){_error_("not implemented yet");};
		void GetInputValue(int* pvalue){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss,IssmDouble time){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss ,int index){_error_("not implemented yet");};
		void GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){_error_("not implemented yet");};
		void GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist);
		IssmDouble InfinityNorm(void);
		IssmDouble Max(void);
		IssmDouble MaxAbs(void);
		IssmDouble Min(void);
		IssmDouble MinAbs(void);
		void Scale(IssmDouble scale_factor);
		void Set(IssmDouble setvalue){_error_("Set not implemented yet");};
		Input* SpawnTriaInput(int index1,int index2,int index3);
		Input* SpawnSegInput(int index1,int index2);
		void SquareMin(IssmDouble* psquaremin,Parameters* parameters);
		void VerticallyIntegrate(Input* thickness_input);

};
#endif  /* _PENTAINPUT_H */
