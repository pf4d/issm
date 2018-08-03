/*!\file: hydrology_core.cpp
 * \brief: core of the hydrology solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void hydrology_core(FemModel* femmodel){

	/*intermediary*/
	int  hydrology_model;
	bool save_results;
	bool modify_loads=true;
	bool isefficientlayer;

	/*first recover parameters common to all solutions*/
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&hydrology_model,HydrologyModelEnum);
	
	if(VerboseSolution()) _printf0_("   computing water heads\n");
			
	/*first compute slopes: */
	if (hydrology_model==HydrologyshreveEnum){
		surfaceslope_core(femmodel);
		bedslope_core(femmodel);
	}

	/*Using the Shreve based Model*/
	if (hydrology_model==HydrologyshreveEnum){
		if(VerboseSolution()) _printf0_("   computing water column\n");
		femmodel->SetCurrentConfiguration(HydrologyShreveAnalysisEnum);
		solutionsequence_nonlinear(femmodel,modify_loads);
		
		/*transfer water column thickness to old water column thickness: */
		InputDuplicatex(femmodel,WatercolumnEnum,WaterColumnOldEnum);
		
		if(save_results){
			if(VerboseSolution()) _printf0_("   saving results \n");
			int outputs[3] = {WatercolumnEnum,HydrologyWaterVxEnum,HydrologyWaterVyEnum};
			femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],3);
			
			/*unload results*/
			if(VerboseSolution()) _printf0_("   saving temporary results\n");
			OutputResultsx(femmodel);
		}
	}

	/*Using the double continuum model*/
	else if (hydrology_model==HydrologydcEnum){
		InputDuplicatex(femmodel,SedimentHeadEnum,SedimentHeadOldEnum);
		femmodel->parameters->FindParam(&isefficientlayer,HydrologydcIsefficientlayerEnum);
		if (isefficientlayer){
			InputDuplicatex(femmodel,EplHeadEnum,EplHeadOldEnum);
			InputDuplicatex(femmodel,HydrologydcEplThicknessEnum,HydrologydcEplThicknessOldEnum);
		}
		
		/*Proceed now to heads computations*/
		solutionsequence_hydro_nonlinear(femmodel);

		if(save_results){
			if(VerboseSolution()) _printf0_("   saving results \n");
			if(isefficientlayer){
				int outputs[9] = {SedimentHeadEnum,SedimentHeadResidualEnum,EplHeadEnum,HydrologydcMaskEplactiveNodeEnum,HydrologydcMaskEplactiveEltEnum,EplHeadSlopeXEnum,EplHeadSlopeYEnum,HydrologydcEplThicknessEnum,EffectivePressureEnum};
				femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],9);
			}
			else{
				int outputs[3] = {SedimentHeadEnum,SedimentHeadResidualEnum,EffectivePressureEnum};
				femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],3);
			}
			/*unload results*/
			if(VerboseSolution()) _printf0_("   saving temporary results\n");
			OutputResultsx(femmodel);
		}
	}

	else if (hydrology_model==HydrologysommersEnum){	
		femmodel->SetCurrentConfiguration(HydrologySommersAnalysisEnum);	
      InputDuplicatex(femmodel,HydrologyHeadEnum,HydrologyHeadOldEnum);	
		solutionsequence_nonlinear(femmodel,modify_loads); 
		if(VerboseSolution()) _printf0_("   updating gap height\n");
		HydrologySommersAnalysis* analysis = new HydrologySommersAnalysis();
		analysis->UpdateGapHeight(femmodel);
		delete analysis;
		
		if(save_results){
			if(VerboseSolution()) _printf0_("   saving results \n");
			int outputs[4] = {HydrologyHeadEnum,HydrologyGapHeightEnum,EffectivePressureEnum,HydrologyBasalFluxEnum};
			femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],4);
			
			/*unload results*/
			if(VerboseSolution()) _printf0_("   saving temporary results\n");
			OutputResultsx(femmodel);
		}
	}

	else{
		_error_("Hydrology model "<< EnumToStringx(hydrology_model) <<" not supported yet");
	}
}

