/*! \file SegInput.h 
 *  \brief: header file for SegInput object
 */

#ifndef _SEGINPUT_H_
#define _SEGINPUT_H_

/*Headers:*/
/*{{{*/
#include "./Input.h"
#include "../Elements/SegRef.h"
class GaussSeg;
class Gauss;
/*}}}*/

class SegInput: public Input,public SegRef{

	public:
		int         enum_type;
		int         interpolation_type;
		IssmDouble* values;

		/*SegInput constructors, destructors*/
		SegInput();
		SegInput(int enum_type,IssmDouble* values,int element_type_in);
		~SegInput();

		/*Object virtual functions definitions*/
		Object *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int     ObjectEnum();

		/*SegInput management:*/
		void   AddTimeValues(IssmDouble* values,int step,IssmDouble time){_error_("not supported yet");};
		void   Configure(Parameters* parameters);
		int  GetResultArraySize(void){_error_("not implemented");};
		int  GetResultInterpolation(void){_error_("not implemented");};
		int  GetResultNumberOfNodes(void){_error_("not implemented");};
		int    InstanceEnum();
		Input* PointwiseDivide(Input* inputB){_error_("not supported yet");};
		Input* PointwiseMax(Input* inputB){_error_("not supported yet");};
		Input* PointwiseMin(Input* inputB){_error_("not supported yet");};
		void ResultToPatch(IssmDouble* values,int nodesperelement,int sid){_error_("not supported yet");};
		Input* SpawnSegInput(int index1,int index2){_error_("not implemented yet");};
		Input* SpawnTriaInput(int index1,int index2,int index3){_error_("not supported yet");};

		/*numerics*/
		void AXPY(Input* xinput,IssmDouble scalar){_error_("not implemented yet");};
		void ChangeEnum(int newenumtype){_error_("not implemented yet");};
		void Constrain(IssmDouble cm_min, IssmDouble cm_max){_error_("not implemented yet");};
		void ConstrainMin(IssmDouble minimum){_error_("not implemented yet");};
		void Extrude(int start){_error_("not supported yet");};
		void GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes){_error_("not implemented yet");};
		void GetInputAverage(IssmDouble* pvalue);
		void GetInputDerivativeAverageValue(IssmDouble* derivativevalues, IssmDouble* xyz_list){_error_("not implemented yet");};
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list,Gauss* gauss);
		void GetInputValue(bool* pvalue){_error_("not implemented yet");}
		void GetInputValue(int* pvalue){_error_("not implemented yet");}
		void GetInputValue(IssmDouble* pvalue){_error_("not implemented yet");}
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss,IssmDouble time){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss ,int index){_error_("not implemented yet");};
		void GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){_error_("not implemented yet");};
		void GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){_error_("not implemented yet");};
		IssmDouble InfinityNorm(void){_error_("not implemented yet");};
		IssmDouble Max(void){_error_("not implemented yet");};
		IssmDouble MaxAbs(void){_error_("not implemented yet");};
		IssmDouble Min(void);
		IssmDouble MinAbs(void){_error_("not implemented yet");};
		void Scale(IssmDouble scale_factor){_error_("not implemented yet");};
		void Set(IssmDouble setvalue){_error_("not implemented yet");};
		void SquareMin(IssmDouble* psquaremin,Parameters* parameters){_error_("not implemented yet");};
		void VerticallyIntegrate(Input* thickness_input){_error_("not supported yet");};

};
#endif  /* _SEGINPUT_H */
