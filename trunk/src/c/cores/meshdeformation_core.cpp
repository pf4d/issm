/*!\file: meshdeformation_core.cpp
 * \brief: core of the slope solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../solutionsequences/solutionsequences.h"
#include "../modules/modules.h"

void meshdeformation_core(FemModel* femmodel){

	/*parameters: */
	bool save_results;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);

	if(VerboseSolution()) _printf0_("computing mesh deformation (elasticity)...\n");

	/*Call on core computations: */
	femmodel->SetCurrentConfiguration(MeshdeformationAnalysisEnum);
	solutionsequence_linear(femmodel);

	if(save_results){
		if(VerboseSolution()) _printf0_("saving results:\n");
	}

}
