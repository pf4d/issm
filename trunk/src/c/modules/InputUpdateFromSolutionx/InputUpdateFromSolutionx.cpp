/*!\file InputUpdateFromSolutionx
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromSolutionx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void InputUpdateFromSolutionx(FemModel* femmodel,Vector<IssmDouble>* solution){

	/*Serialize solution, so that elements can index into it on every CPU: */
	IssmDouble* serial_solution=solution->ToMPISerial();

	/*Call overloaded form of InputUpdateFromSolutionx: */
	InputUpdateFromSolutionx(femmodel,serial_solution);

	/*cleanup and return*/
	xDelete<IssmDouble>(serial_solution);
}

void InputUpdateFromSolutionx(FemModel* femmodel,IssmDouble* solution){

	/*retrive parameters: */
	int analysisenum;
	femmodel->parameters->FindParam(&analysisenum,AnalysisTypeEnum);

	/*Display message*/
	if(VerboseModule()) _printf0_("   Updating inputs from solution\n");

	Analysis* analysis = EnumToAnalysis(analysisenum);
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		analysis->InputUpdateFromSolution(solution,element);
	}
	delete analysis;
}
