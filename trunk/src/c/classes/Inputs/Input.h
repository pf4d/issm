/*!\file:  Input.h
 * \brief abstract class for Input object
 */ 

#ifndef _INPUT_H_
#define _INPUT_H_

/*Headers:*/
/*{{{*/
#include "../../datastructures/datastructures.h"
#include "../../shared/shared.h"
class Node;
class Gauss;
class Gauss;
class GaussSeg;
class Parameters;
class Gauss;
template <class doubletype> class Vector;
/*}}}*/

class Input: public Object{

	public: 

		virtual        ~Input(){};

		virtual void ChangeEnum(int newenumtype)=0;
		virtual void Configure(Parameters* parameters)=0;
		virtual void GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes)=0;
		virtual void GetInputAverage(IssmDouble* pvalue)=0;
		virtual void GetInputDerivativeAverageValue(IssmDouble* derivativevalues, IssmDouble* xyz_list)=0;
		virtual void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list, Gauss* gauss)=0;
		virtual void GetInputValue(bool* pvalue)=0;
		virtual void GetInputValue(int* pvalue)=0;
		virtual void GetInputValue(IssmDouble* pvalue)=0;
		virtual void GetInputValue(IssmDouble* pvalue,Gauss* gauss)=0;
		virtual void GetInputValue(IssmDouble* pvalue,Gauss* gauss,IssmDouble time)=0;
		virtual void GetInputValue(IssmDouble* pvalue,Gauss* gauss,int index)=0;
		virtual void GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime)=0;
		virtual int  InstanceEnum()=0; 

		virtual void   AXPY(Input* xinput,IssmDouble scalar)=0;
		virtual void   Constrain(IssmDouble cm_min, IssmDouble cm_max)=0;
		virtual void   ConstrainMin(IssmDouble minimum)=0;
		virtual void   Extrude(int start)=0;
		virtual void   GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist)=0;
		virtual IssmDouble InfinityNorm(void)=0;
		virtual IssmDouble Max(void)=0;
		virtual IssmDouble MaxAbs(void)=0;
		virtual IssmDouble Min(void)=0;
		virtual IssmDouble MinAbs(void)=0;
		virtual void   Scale(IssmDouble scale_factor)=0;
		virtual void   Set(IssmDouble setvalue)=0;
		virtual void   SquareMin(IssmDouble* psquaremin,Parameters* parameters)=0;
		virtual void   VerticallyIntegrate(Input* thickness_input)=0;

		virtual int  GetResultArraySize(void)=0;
		virtual int  GetResultInterpolation(void)=0;
		virtual int  GetResultNumberOfNodes(void)=0;
		virtual Input* PointwiseDivide(Input* inputB)=0;
		virtual Input* PointwiseMax(Input* inputmax)=0;
		virtual Input* PointwiseMin(Input* inputmin)=0;
		virtual Input* SpawnSegInput(int index1,int index2)=0;
		virtual Input* SpawnTriaInput(int index1,int index2,int index3)=0;
		virtual void ResultToMatrix(IssmDouble* values,int ncols,int sid){_error_("not supported yet");};
		virtual void ResultToPatch(IssmDouble* values,int nodesperelement,int sid){_error_("not supported yet");}; 
};
#endif
