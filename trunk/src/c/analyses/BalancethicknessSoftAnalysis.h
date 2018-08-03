/*! \file BalancethicknessSoftAnalysis.h 
 *  \brief: header file for generic external result object
 */

#ifndef _BalancethicknessSoftAnalysis_
#define _BalancethicknessSoftAnalysis_

/*Headers*/
#include "./Analysis.h"

class BalancethicknessSoftAnalysis: public Analysis{

	public:
		/*Model processing*/
		int  DofsPerNode(int** doflist,int domaintype,int approximation);
		void UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum);
		void UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type);
		void CreateNodes(Nodes* nodes,IoModel* iomodel);
		void CreateConstraints(Constraints* constraints,IoModel* iomodel);
		void CreateLoads(Loads* loads, IoModel* iomodel);

		/*Finite element Analysis*/
		void           Core(FemModel* femmodel);
		ElementVector* CreateDVector(Element* element);
		ElementMatrix* CreateJacobianMatrix(Element* element);
		ElementMatrix* CreateKMatrix(Element* element);
		ElementVector* CreatePVector(Element* element);
		void GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element);
		void GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index);
		void InputUpdateFromSolution(IssmDouble* solution,Element* element);
		void UpdateConstraints(FemModel* femmodel);
};
#endif
