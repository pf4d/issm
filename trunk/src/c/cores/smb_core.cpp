/*!\file: smb_core.cpp
 * \brief: core of the smb solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void smb_core(FemModel* femmodel){

	/*parameters: */
	Analysis* analysis=NULL;
	int    smb_model;
	int    numoutputs;
	bool   save_results;
	int    solution_type;
	char** requested_outputs = NULL;

	/*activate configuration*/
	femmodel->SetCurrentConfiguration(SmbAnalysisEnum);

	/*recover parameters: */
	femmodel->parameters->FindParam(&smb_model,SmbEnum);
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&numoutputs,SmbNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,SmbRequestedOutputsEnum);
			
	if(VerboseSolution()) _printf0_("   computing smb \n");
 
	analysis = new SmbAnalysis();
	analysis->Core(femmodel);
	delete analysis;

	if(save_results){
		if(VerboseSolution()) _printf0_("   saving results\n");
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}

	if(solution_type==SmbSolutionEnum)femmodel->RequestedDependentsx();

	/*Free ressources:*/
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}
}
