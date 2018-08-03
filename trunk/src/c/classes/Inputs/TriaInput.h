/*! \file TriaInput.h 
 *  \brief: header file for TriaInput object
 */

#ifndef _TRIAINPUT_H_
#define _TRIAINPUT_H_

/*Headers:*/
/*{{{*/
#include "./Input.h"
#include "../Elements/TriaRef.h"
class Gauss;
class Gauss;
/*}}}*/

class TriaInput: public Input,public TriaRef{

	public:
		int         enum_type;
		int         interpolation_type;
		IssmDouble* values;

		/*TriaInput constructors, destructors*/
		TriaInput();
		TriaInput(int enum_type,IssmDouble* values,int element_type_in);
		~TriaInput();

		/*Object virtual functions definitions*/
		Object *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int     ObjectEnum();

		/*TriaInput management:*/
		void   AddTimeValues(IssmDouble* values,int step,IssmDouble time){_error_("not supported yet");};
		void   Configure(Parameters* parameters);
		int    GetResultArraySize(void);
		int    GetResultInterpolation(void);
		int    GetResultNumberOfNodes(void);
		int    InstanceEnum();
		Input* PointwiseDivide(Input* inputB);
		Input* PointwiseMax(Input* inputB);
		Input* PointwiseMin(Input* inputB);
		void   ResultToPatch(IssmDouble* values,int nodesperelement,int sid);
		Input* SpawnSegInput(int index1,int index2);
		Input* SpawnTriaInput(int index1,int index2,int index3);

		/*numerics*/
		void AXPY(Input* xinput,IssmDouble scalar);
		void ChangeEnum(int newenumtype);
		void Constrain(IssmDouble cm_min, IssmDouble cm_max);
		void ConstrainMin(IssmDouble minimum);
		void Extrude(int start){_error_("not supported yet");};
		void GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes);
		void GetInputAverage(IssmDouble* pvalue);
		void GetInputDerivativeAverageValue(IssmDouble* derivativevalues, IssmDouble* xyz_list);
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list,Gauss* gauss);
		void GetInputValue(bool* pvalue){_error_("not implemented yet");}
		void GetInputValue(int* pvalue){_error_("not implemented yet");}
		void GetInputValue(IssmDouble* pvalue){_error_("not implemented yet (Input is "<<EnumToStringx(this->enum_type)<<")");}//{_error_("not implemented yet");}
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss,IssmDouble time){_error_("not implemented yet");};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss,int index){_error_("not implemented yet");};
		void GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime);
		void GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist);
		IssmDouble InfinityNorm(void);
		IssmDouble Max(void);
		IssmDouble MaxAbs(void);
		IssmDouble Min(void);
		IssmDouble MinAbs(void);
		void Scale(IssmDouble scale_factor);
		void Set(IssmDouble setvalue);
		void SquareMin(IssmDouble* psquaremin,Parameters* parameters);
		void VerticallyIntegrate(Input* thickness_input){_error_("not supported yet");};

};
#endif  /* _TRIAINPUT_H */
