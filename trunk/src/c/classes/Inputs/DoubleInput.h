/*! \file DoubleInput.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _DOUBLEINPUT_H_
#define _DOUBLEINPUT_H_

/*Headers:*/
/*{{{*/
#include "./Input.h"
class Gauss;
class Gauss;
/*}}}*/

class DoubleInput: public Input{

	public:
		int    enum_type;
		IssmDouble value;

		/*DoubleInput constructors, destructors: {{{*/
		DoubleInput();
		DoubleInput(int enum_type,IssmDouble value);
		~DoubleInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int   ObjectEnum();
		/*}}}*/
		/*DoubleInput management: {{{*/
		void AddTimeValues(IssmDouble* values,int step,IssmDouble time){_error_("not supported yet");};
		void Configure(Parameters* parameters);
		int  GetResultArraySize(void){return 1;};
		int  GetResultInterpolation(void){return P0Enum;};
		int  GetResultNumberOfNodes(void){return 1;};
		int   InstanceEnum();
		Input* PointwiseDivide(Input* inputB);
		Input* PointwiseMax(Input* inputB);
		Input* PointwiseMin(Input* inputB);
		void ResultToPatch(IssmDouble* values,int nodesperelement,int sid){_error_("not supported yet");};
		Input* SpawnSegInput(int index1,int index2);
		Input* SpawnTriaInput(int index1,int index2,int index3);
		/*}}}*/
		/*numerics: {{{*/
		void AXPY(Input* xinput,IssmDouble scalar);
		void ChangeEnum(int newenumtype);
		void Constrain(IssmDouble cm_min, IssmDouble cm_max);
		void ConstrainMin(IssmDouble minimum);
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
		IssmDouble InfinityNorm(void){_error_("not implemented yet");};
		IssmDouble Max(void);
		IssmDouble MaxAbs(void);
		IssmDouble Min(void);
		IssmDouble MinAbs(void);
		void Scale(IssmDouble scale_factor);
		void Set(IssmDouble setvalue){_error_("Set not implemented yet");};
		void SquareMin(IssmDouble* psquaremin,Parameters* parameters);
		void VerticallyIntegrate(Input* thickness_input);
		/*}}}*/

};
#endif  /* _DOUBLEINPUT_H */
