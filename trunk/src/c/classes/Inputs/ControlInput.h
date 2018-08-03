/*! \file ControlInput.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _CONTROLINPUT_H_
#define _CONTROLINPUT_H_

/*Headers:*/
/*{{{*/
#include "./Input.h"
class Gauss;
class Gauss;
/*}}}*/

class ControlInput: public Input{

	public:
		int    control_id;
		int    enum_type;
		Input* gradient;
		Input* maxvalues;
		Input* minvalues;
		Input* savedvalues;
		Input* values;

		/*ControlInput constructors, destructors: {{{*/
		ControlInput();
		ControlInput(int enum_type,int enum_input,IssmDouble* pvalues,IssmDouble* pmin,IssmDouble* pmax,int id);
		~ControlInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int   ObjectEnum();
		/*}}}*/
		/*ControlInput management: {{{*/
		void AddTimeValues(IssmDouble* values,int step,IssmDouble time){_error_("not supported yet");};
		void Configure(Parameters* parameters);
		Input* PointwiseDivide(Input* inputB){_error_("not implemented yet");};
		Input* PointwiseMax(Input* inputB){_error_("not implemented yet");};
		Input* PointwiseMin(Input* inputB){_error_("not implemented yet");};
		int    InstanceEnum();
		Input* SpawnSegInput(int index1,int index2);
		Input* SpawnTriaInput(int index1,int index2,int index3);
		/*}}}*/
		/*numerics: {{{*/
		void AXPY(Input* xinput,IssmDouble scalar);
		void Constrain(void);
		void Constrain(IssmDouble min,IssmDouble max);
		void ConstrainMin(IssmDouble minimum){_error_("not implemented yet");};
		void ChangeEnum(int newenumtype){_error_("not implemented yet");};
		void Extrude(int start);
		void GetGradient(Vector<IssmDouble>* gradient_vec,int* doflist);
		void GetGradientValue(IssmDouble* pvalue,Gauss* gauss);
		void GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes){_error_("not implemented yet");};
		void GetInputAverage(IssmDouble* pvalue);
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list,Gauss* gauss);
		void GetInputDerivativeAverageValue(IssmDouble* derivativevalues, IssmDouble* xyz_list){_error_("not implemented yet");};
		void GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){_error_("not implemented yet");};
		void GetInputValue(bool* pvalue);
		void GetInputValue(int* pvalue);
		void GetInputValue(IssmDouble* pvalue);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss,IssmDouble time){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss ,int index){_error_("not implemented yet");};
		int  GetResultArraySize(void){return 1;};
		int  GetResultInterpolation(void);
		int  GetResultNumberOfNodes(void);
		void GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist,const char* data);
		void GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist);
		IssmDouble InfinityNorm(void){_error_("not implemented yet");};
		IssmDouble Max(void){_error_("not implemented yet");};
		IssmDouble MaxAbs(void){_error_("not implemented yet");};
		IssmDouble Min(void);
		IssmDouble MinAbs(void){_error_("not implemented yet");};
		void ResultToPatch(IssmDouble* values,int nodesperelement,int sid){_error_("not supported yet");};
		void SaveValue(void);
		void Scale(IssmDouble scale_factor){_error_("not implemented yet");};
		void ScaleGradient(IssmDouble scale);
		void Set(IssmDouble setvalue){_error_("Set not implemented yet");};
		void SetGradient(Input* gradient_in);
		void SetInput(Input* in_input);
		void SquareMin(IssmDouble* psquaremin,Parameters* parameters){_error_("not implemented yet");};
		void UpdateValue(IssmDouble scalar);
		void VerticallyIntegrate(Input* thickness_input);
		/*}}}*/

};
#endif  /* _CONTROLINPUT_H */
