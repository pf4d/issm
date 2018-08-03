/*! \file IntInput.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _INTINPUT_H_
#define _INTINPUT_H_

/*Headers:*/
/*{{{*/
#include "./Input.h"
class Gauss;
class Gauss;
/*}}}*/

class IntInput: public Input{

	public:
		/*just hold 3 values for 3 vertices: */
		int    enum_type;
		IssmInt value;

		/*IntInput constructors, destructors: {{{*/
		IntInput();
		IntInput(int enum_type,IssmInt value);
		~IntInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int   ObjectEnum();
		/*}}}*/
		/*IntInput management: {{{*/
		void AddTimeValues(IssmDouble* values,int step,IssmDouble time){_error_("not supported yet");};
		void Configure(Parameters* parameters);
		int  GetResultArraySize(void){return 1;};
		int  GetResultInterpolation(void){return P0Enum;};
		int  GetResultNumberOfNodes(void){return 1;};
		int   InstanceEnum();
		Input* PointwiseDivide(Input* inputB){_error_("not implemented yet");};
		Input* PointwiseMax(Input* inputB){_error_("not implemented yet");};
		Input* PointwiseMin(Input* inputB){_error_("not implemented yet");};
		void ResultToPatch(IssmDouble* values,int nodesperelement,int sid){_error_("not supported yet");};
		Input* SpawnSegInput(int index1,int index2);
		Input* SpawnTriaInput(int index1,int index2,int index3);
		/*}}}*/
		/*numerics: {{{*/
		void AXPY(Input* xinput,IssmDouble scalar);
		void ChangeEnum(int newenumtype);
		void Constrain(IssmDouble cm_min, IssmDouble cm_max);
		void ConstrainMin(IssmDouble minimum){_error_("not implemented yet");};
		void Extrude(int start){_error_("not supported yet");};
		void GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes){_error_("not implemented yet");};
		void GetInputAverage(IssmDouble* pvalue);
		void GetInputDerivativeAverageValue(IssmDouble* derivativevalues, IssmDouble* xyz_list){_error_("not implemented yet");};
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void GetInputValue(bool* pvalue);
		void GetInputValue(int* pvalue);
		void GetInputValue(IssmDouble* pvalue);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss,IssmDouble time){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss ,int index){_error_("not implemented yet");};
		void GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){_error_("not implemented yet");};
		void GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist);
		IssmDouble InfinityNorm(void){_error_("InfinityNorm not implemented for integers");};
		IssmDouble Max(void){_error_("Max not implemented for integers");};
		IssmDouble MaxAbs(void){_error_("Max not implemented for integers");};
		IssmDouble Min(void){_error_("Min not implemented for integers");};
		IssmDouble MinAbs(void){_error_("Min not implemented for integers");};
		void Scale(IssmDouble scale_factor);
		void Set(IssmDouble setvalue){_error_("Set not implemented yet");};
		void SquareMin(IssmDouble* psquaremin,Parameters* parameters);
		void VerticallyIntegrate(Input* thickness_input){_error_("not supported yet");};
		/*}}}*/

};
#endif  /* _INTINPUT_H */
