/*! \file DoubleArrayInput.h 
 *  \brief: header file for vector type input object
 */

#ifndef _DOUBLE_ARRAY_INPUT_H_
#define _DOUBLE_ARRAY_INPUT_H_

/*Headers:*/
/*{{{*/
#include "./Input.h"
/*}}}*/

class DoubleArrayInput: public Input{

	public:
		int    enum_type;
		IssmDouble* values; /*vector*/
		int         m; /*size of vector*/

		/*DoubleArrayInput constructors, destructors: {{{*/
		DoubleArrayInput();
		DoubleArrayInput(int enum_type,IssmDouble* values, int m);
		~DoubleArrayInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int   ObjectEnum();
		/*}}}*/
		/*DoubleArrayInput management: {{{*/
		void Configure(Parameters* parameters);
		void GetValues(IssmDouble** pvalues,int* pm);
		int  GetResultArraySize(void){return m;};
		int  GetResultInterpolation(void){return P0ArrayEnum;};
		int  GetResultNumberOfNodes(void){return 1;};
		int   InstanceEnum();
		Input* PointwiseDivide(Input* inputB){_error_("not implemented yet");};
		Input* PointwiseMax(Input* inputB){_error_("not implemented yet");};
		Input* PointwiseMin(Input* inputB){_error_("not implemented yet");};
		void ResultToMatrix(IssmDouble* values,int ncols,int sid);
		Input* SpawnSegInput(int index1,int index2){_error_("not implemented yet");};
		Input* SpawnTriaInput(int index1,int index2,int index3){_error_("not implemented yet");};
		/*}}}*/
		/*numerics: {{{*/
		void AXPY(Input* xinput,IssmDouble scalar){_error_("not implemented yet");};
		void ChangeEnum(int newenumtype);
		void Constrain(IssmDouble cm_min, IssmDouble cm_max){_error_("not implemented yet");};
		void ConstrainMin(IssmDouble minimum){_error_("not implemented yet");};
		void Extrude(int start){_error_("not supported yet");};
		void GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes){_error_("not implemented yet");};
		void GetInputAverage(IssmDouble* pvalue){_error_("not implemented yet");};
		void GetInputDerivativeAverageValue(IssmDouble* derivativevalues, IssmDouble* xyz_list){_error_("not implemented yet");};
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void GetInputValue(bool* pvalue){_error_("not implemented yet");};
		void GetInputValue(int* pvalue){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss,IssmDouble time){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss ,int index){_error_("not implemented yet");};
		void GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){_error_("not implemented yet");};
		void GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){_error_("not implemented yet");};
		IssmDouble InfinityNorm(void){_error_("not implemented yet");};
		IssmDouble Max(void){_error_("not implemented yet");};
		IssmDouble MaxAbs(void){_error_("not implemented yet");};
		IssmDouble Min(void){_error_("not implemented yet");};
		IssmDouble MinAbs(void){_error_("not implemented yet");};
		void Scale(IssmDouble scale_factor){_error_("not implemented yet");};
		void Set(IssmDouble setvalue){_error_("Set not implemented yet");};
		void SquareMin(IssmDouble* psquaremin,Parameters* parameters){_error_("not implemented yet");};
		void VerticallyIntegrate(Input* thickness_input){_error_("not implemented yet");};
		/*}}}*/

};
#endif  /* _DOUBLE_ARRAY_INPUT_H */
