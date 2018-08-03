/*!\file: balancethickness_core.cpp
 * \brief: core of the balancethickness solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void balancethickness2_core(FemModel* femmodel){

	/*parameters: */
	bool        save_results;
	IssmDouble  l = 3.;

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);

	if(VerboseSolution()) _printf0_("computing smooth surface slopes:\n");
	//femmodel->parameters->SetParam(l,SmoothThicknessMultiplierEnum);
	//femmodel->SetCurrentConfiguration(SmoothAnalysisEnum);
	//femmodel->parameters->SetParam(SurfaceSlopeXEnum,InputToSmoothEnum);
	//solutionsequence_linear(femmodel);
	//femmodel->parameters->SetParam(SurfaceSlopeYEnum,InputToSmoothEnum);
	//solutionsequence_linear(femmodel);
	surfaceslope_core(femmodel);

	if(VerboseSolution()) _printf0_("call computational core:\n");
	femmodel->SetCurrentConfiguration(Balancethickness2AnalysisEnum);
	solutionsequence_linear(femmodel);
	//solutionsequence_nonlinear(femmodel,false);

	if(save_results){
		if(VerboseSolution()) _printf0_("   saving results\n");
		const int numoutputs = 6;
		int outputs[numoutputs] = {SurfaceEnum,SurfaceSlopeXEnum,SurfaceSlopeYEnum,VxEnum,VyEnum,VelEnum};
		//const int numoutputs = 4;
		//int outputs[numoutputs] = {SurfaceEnum,VxEnum,VyEnum,VelEnum};
		femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],numoutputs);
	}

}
