/*!\file: balancethickness_core.cpp
 * \brief: core of the balancethickness solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void balancethickness_core(FemModel* femmodel){

	/*parameters: */
	bool save_results;

	/*activate formulation: */
	femmodel->SetCurrentConfiguration(BalancethicknessAnalysisEnum);

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);

	if(VerboseSolution()) _printf0_("call computational core:\n");
	solutionsequence_linear(femmodel);

	if(save_results){
		if(VerboseSolution()) _printf0_("   saving results\n");
		int outputs = ThicknessEnum;
		femmodel->RequestedOutputsx(&femmodel->results,&outputs,1);
	}

}
