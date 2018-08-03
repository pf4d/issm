/*!\file: masstransport_core.cpp
 * \brief: core of the masstransport solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void masstransport_core(FemModel* femmodel){

	/*parameters: */
	int    numoutputs,domaintype;
	bool   save_results;
	bool   isFS,isfreesurface,dakota_analysis;
	int    solution_type,stabilization;
	char** requested_outputs = NULL;

	/*activate configuration*/
	femmodel->SetCurrentConfiguration(MasstransportAnalysisEnum);

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&isFS,FlowequationIsFSEnum);
	femmodel->parameters->FindParam(&isfreesurface,MasstransportIsfreesurfaceEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&numoutputs,MasstransportNumRequestedOutputsEnum);
	femmodel->parameters->FindParam(&stabilization,MasstransportStabilizationEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,MasstransportRequestedOutputsEnum);
			
	if(VerboseSolution()) _printf0_("   computing mass transport\n");

	/*Transport mass or free surface*/
	if(isFS && isfreesurface){
		if(VerboseSolution()) _printf0_("   call free surface computational core\n");
		femmodel->SetCurrentConfiguration(FreeSurfaceBaseAnalysisEnum);
		solutionsequence_linear(femmodel);
		femmodel->parameters->SetParam(BaseEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);

		femmodel->SetCurrentConfiguration(FreeSurfaceTopAnalysisEnum);
		solutionsequence_linear(femmodel);
		femmodel->parameters->SetParam(SurfaceEnum,InputToExtrudeEnum);
		extrudefromtop_core(femmodel);
	}
	else{
		if(VerboseSolution()) _printf0_("   call computational core\n");
		femmodel->parameters->SetParam(VxEnum,InputToDepthaverageInEnum);
		femmodel->parameters->SetParam(VxAverageEnum,InputToDepthaverageOutEnum);
		depthaverage_core(femmodel);
		if(domaintype!=Domain2DverticalEnum){
			femmodel->parameters->SetParam(VyEnum,InputToDepthaverageInEnum);
			femmodel->parameters->SetParam(VyAverageEnum,InputToDepthaverageOutEnum);
			depthaverage_core(femmodel);
		}
		femmodel->SetCurrentConfiguration(MasstransportAnalysisEnum);
		if(stabilization==4){
			solutionsequence_fct(femmodel);
		}
		else{
			solutionsequence_linear(femmodel);
		}
		femmodel->parameters->SetParam(ThicknessEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);
		femmodel->parameters->SetParam(BaseEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);
		femmodel->parameters->SetParam(SurfaceEnum,InputToExtrudeEnum);
		extrudefrombase_core(femmodel);
	}

	if(save_results){
		if(VerboseSolution()) _printf0_("   saving results\n");
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}

	if(solution_type==MasstransportSolutionEnum)femmodel->RequestedDependentsx();

	/*Free ressources:*/
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}
}
