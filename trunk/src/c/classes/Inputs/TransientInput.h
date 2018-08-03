/*! \file TransientInput.h 
 *  \brief: header file for transientinput object
 */

#ifndef _TRANSIENTINPUT_H_
#define _TRANSIENTINPUT_H_

/*Headers:*/
/*{{{*/
#include "./Input.h"
class Gauss;
class Parameters;
class Gauss;
/*}}}*/

class TransientInput: public Input{

	public:
		int         enum_type;
		int         numtimesteps;
		Inputs     *inputs;
		IssmDouble *timesteps;
		Parameters *parameters;     //to find current time.

		/*TransientInput constructors, destructors: {{{*/
		TransientInput();
		TransientInput(int enum_type);
		TransientInput(int in_enum_type,IssmDouble* times,int N);
		~TransientInput();
		void AddTimeInput(Input* input,IssmDouble time);
		void AddTimeInput(Input* input);
		/*}}}*/
		/*Object virtual functions definitions:{{{*/
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id();
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int   ObjectEnum();
		/*}}}*/
		/*TransientInput management: {{{*/
		void Configure(Parameters* parameters);
		int  GetResultArraySize(void);
		int  GetResultInterpolation(void);
		int  GetResultNumberOfNodes(void);
		int    InstanceEnum();
		Input* PointwiseDivide(Input* input_in){_error_("not implemented yet");};
		Input* PointwiseMax(Input* input_in){_error_("not implemented yet");};
		Input* PointwiseMin(Input* input_in){_error_("not implemented yet");};
		void ResultToPatch(IssmDouble* values,int nodesperelement,int sid){_error_("not supported yet");};
		Input* SpawnSegInput(int index1,int index2);
		Input* SpawnTriaInput(int index1,int index2,int index3);
		/*}}}*/
		/*numerics: {{{*/
		void AXPY(Input* xforcing,IssmDouble scalar){_error_("not implemented yet");};
		void ChangeEnum(int newenumtype);
		void Constrain(IssmDouble cm_min, IssmDouble cm_max){_error_("not implemented yet");};
		void ConstrainMin(IssmDouble minimum){_error_("not implemented yet");};
		void Extrude(int start);
		void GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes);
		void GetInputAverage(IssmDouble* pvalue);
		void GetInputDerivativeAverageValue(IssmDouble* derivativevalues, IssmDouble* xyz_list){_error_("not implemented yet");};
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list,Gauss* gauss);
		void GetInputValue(bool* pvalue){_error_("not implemented yet");};
		void GetInputValue(int* pvalue){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss,IssmDouble time);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss ,int index){_error_("not implemented yet");};
		void GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime);
		Input* GetTimeInput(IssmDouble time);
		void GetTimeValues(IssmDouble* values,IssmDouble time){_error_("not implemented yet");};
		void GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist);
		IssmDouble InfinityNorm(void);
		IssmDouble Max(void);
		IssmDouble MaxAbs(void);
		IssmDouble Min(void);
		IssmDouble MinAbs(void);
		void Scale(IssmDouble scale_factor){_error_("not implemented yet");};
		void Set(IssmDouble setvalue){_error_("Set not implemented yet");};
		void SquareMin(IssmDouble* psquaremin,Parameters* parameters);
		void VerticallyIntegrate(Input* thickness_forcing){_error_("not supported yet");};
		/*}}}*/

};
#endif  /* _TRANSIENTINPUT_H */
