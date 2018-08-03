/*! \file DatasetInput.h 
 *  \brief: header file for datasetinput object
 */

#ifndef _DATASETINPUT_H_
#define _DATASETINPUT_H_

/*Headers:*/
/*{{{*/
#include "./Input.h"
class Gauss;
class Gauss;
/*}}}*/

class DatasetInput: public Input{

	public:
		int     enum_type;
		int     numids;
		Inputs *inputs;
		int    *ids;

		/*DatasetInput constructors, destructors: {{{*/
		DatasetInput();
		DatasetInput(int enum_type);
		~DatasetInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id();
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int   ObjectEnum();
		/*}}}*/
		/*DatasetInput management: {{{*/
		void   AddInput(Input* input,int id);
		void AddTimeValues(IssmDouble* values,int step,IssmDouble time){_error_("not supported yet");};
		void Configure(Parameters* parameters);
		int    InstanceEnum();
		Input* PointwiseDivide(Input* inputB){_error_("not implemented yet");};
		Input* PointwiseMin(Input* inputB){_error_("not implemented yet");};
		Input* PointwiseMax(Input* inputB){_error_("not implemented yet");};
		Input* SpawnSegInput(int index1,int index2);
		Input* SpawnTriaInput(int index1,int index2,int index3);
		/*}}}*/
		/*numerics: {{{*/
		void AXPY(Input* xinput,IssmDouble scalar){_error_("not implemented yet");};
		void ChangeEnum(int newenumtype){_error_("not implemented yet");};
		void Constrain(void){_error_("not implemented yet");};
		void Constrain(IssmDouble min,IssmDouble max){_error_("not implemented yet");};
		void ConstrainMin(IssmDouble minimum){_error_("not implemented yet");};
		void Extrude(int start){_error_("not implemented yet");};
		void GetGradient(Vector<IssmDouble>* gradient_vec,int* doflist){_error_("not implemented yet");};
		void GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes){_error_("not implemented yet");};
		void GetInputAverage(IssmDouble* pvalue){_error_("not implemented yet");};
		void GetInputDerivativeAverageValue(IssmDouble* derivativevalues, IssmDouble* xyz_list){_error_("not implemented yet");};
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void GetInputValue(bool* pvalue){_error_("not implemented yet");};
		void GetInputValue(int* pvalue){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss,IssmDouble time){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss ,int index);
		void GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){_error_("not implemented yet");};
		int GetResultArraySize(void){_error_("not implemented yet");};
		int GetResultInterpolation(void){_error_("not implemented yet");};
		int GetResultNumberOfNodes(void){_error_("not implemented yet");};
		void GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){_error_("not implemented yet");};
		IssmDouble InfinityNorm(void){_error_("not implemented yet");};
		IssmDouble Max(void){_error_("not implemented yet");};
		IssmDouble MaxAbs(void){_error_("not implemented yet");};
		IssmDouble Min(void){_error_("not implemented yet");};
		IssmDouble MinAbs(void){_error_("not implemented yet");};
		void ResultToPatch(IssmDouble* values,int nodesperelement,int sid){_error_("not supported yet");};
		void SaveValue(void){_error_("not implemented yet");};
		void Scale(IssmDouble scale_factor){_error_("not implemented yet");};
		void ScaleGradient(IssmDouble scale){_error_("not implemented yet");};
		void Set(IssmDouble setvalue){_error_("Set not implemented yet");};
		void SetGradient(Input* gradient_in){_error_("not implemented yet");};
		void SquareMin(IssmDouble* psquaremin,Parameters* parameters){_error_("not implemented yet");};
		void UpdateValue(IssmDouble scalar){_error_("not implemented yet");};
		void VerticallyIntegrate(Input* thickness_input){_error_("not implemented yet");};
		/*}}}*/

};
#endif  /* _DATASETINPUT_H */
