/*!\file: transient_3d_core.cpp
 * \brief: core of the transient_3d solution 
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <float.h>
#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void transient_core(FemModel* femmodel){

	/*parameters: */
	IssmDouble finaltime,dt,yts;
	bool       isstressbalance,ismasstransport,issmb,isFS,isthermal,isgroundingline,isgia,isslr,iscoupler,ismovingfront,isdamageevolution,ishydrology;
	bool       save_results,dakota_analysis;
	bool       time_adapt;
	int        output_frequency;
	int        recording_frequency;
	int        domaintype,groundingline_migration,smb_model,amr_frequency;
	int        numoutputs;
	Analysis  *analysis          = NULL;
	char     **requested_outputs = NULL;

	/*intermediary: */
	int        step;
	IssmDouble time;

	/*first, figure out if there was a check point, if so, do a reset of the FemModel* femmodel structure. */
	femmodel->parameters->FindParam(&recording_frequency,SettingsRecordingFrequencyEnum);
	if(recording_frequency) femmodel->Restart();

	/*then recover parameters common to all solutions*/
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&step,StepEnum);
	femmodel->parameters->FindParam(&time,TimeEnum);
	femmodel->parameters->FindParam(&finaltime,TimesteppingFinalTimeEnum);
	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);
	femmodel->parameters->FindParam(&dakota_analysis,QmuIsdakotaEnum);
	femmodel->parameters->FindParam(&output_frequency,SettingsOutputFrequencyEnum);
	femmodel->parameters->FindParam(&time_adapt,TimesteppingTimeAdaptEnum);
	femmodel->parameters->FindParam(&isstressbalance,TransientIsstressbalanceEnum);
	femmodel->parameters->FindParam(&ismasstransport,TransientIsmasstransportEnum);
	femmodel->parameters->FindParam(&issmb,TransientIssmbEnum);
	femmodel->parameters->FindParam(&isthermal,TransientIsthermalEnum);
	femmodel->parameters->FindParam(&isgia,TransientIsgiaEnum);
	femmodel->parameters->FindParam(&isslr,TransientIsslrEnum);
	femmodel->parameters->FindParam(&iscoupler,TransientIscouplerEnum);
	femmodel->parameters->FindParam(&isgroundingline,TransientIsgroundinglineEnum);
	femmodel->parameters->FindParam(&ismovingfront,TransientIsmovingfrontEnum);
	femmodel->parameters->FindParam(&isdamageevolution,TransientIsdamageevolutionEnum);
	femmodel->parameters->FindParam(&ishydrology,TransientIshydrologyEnum);
	femmodel->parameters->FindParam(&amr_frequency,TransientAmrFrequencyEnum);
	femmodel->parameters->FindParam(&isFS,FlowequationIsFSEnum);
	if(isgroundingline) femmodel->parameters->FindParam(&groundingline_migration,GroundinglineMigrationEnum);
	femmodel->parameters->FindParam(&numoutputs,TransientNumRequestedOutputsEnum);
	if(numoutputs) femmodel->parameters->FindParam(&requested_outputs,&numoutputs,TransientRequestedOutputsEnum);

	#ifdef _HAVE_NEOPZ_
	bool ismismip = false;//itapopo testing restart 
	if(ismismip) femmodel->ReMesh();
	#endif

	while(time < finaltime - (yts*DBL_EPSILON)){ //make sure we run up to finaltime.
		/*Increment*/
		if(time_adapt){
			femmodel->TimeAdaptx(&dt);
			if(time+dt>finaltime) dt=finaltime-time;
			femmodel->parameters->SetParam(dt,TimesteppingTimeStepEnum);
		}
		step+=1;
		time+=dt;
		femmodel->parameters->SetParam(time,TimeEnum);
		femmodel->parameters->SetParam(step,StepEnum);

		if(VerboseSolution()) _printf0_("iteration " << step << "/" << floor((finaltime-time)/dt)+step << "  time [yr]: " << setprecision(4) << time/yts << " (time step: " << dt/yts << ")\n");
		if(step%output_frequency==0 || (time >= finaltime - (yts*DBL_EPSILON)) || step==1)
		 save_results=true;
		else
		 save_results=false;
		femmodel->parameters->SetParam(save_results,SaveResultsEnum);

		if(isthermal && domaintype==Domain3DEnum){ 
			if(issmb){
				bool isenthalpy;
				femmodel->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);
				femmodel->parameters->FindParam(&smb_model,SmbEnum);
				if(isenthalpy){
					if(smb_model==SMBpddEnum)     ResetBoundaryConditions(femmodel,EnthalpyAnalysisEnum);
					if(smb_model==SMBd18opddEnum) ResetBoundaryConditions(femmodel,EnthalpyAnalysisEnum);
				}
				else{
					if(smb_model==SMBpddEnum)     ResetBoundaryConditions(femmodel,ThermalAnalysisEnum);
					if(smb_model==SMBd18opddEnum) ResetBoundaryConditions(femmodel,ThermalAnalysisEnum);
				}
			}
			if(VerboseSolution()) _printf0_("   computing thermal regime\n");
			thermal_core(femmodel);
		}

		if(ishydrology) hydrology_core(femmodel);

		if(isstressbalance) stressbalance_core(femmodel);

		if(isdamageevolution) damage_core(femmodel);

		if(ismovingfront)	movingfront_core(femmodel);

		/* from here on, prepare geometry for next time step*/
		if(issmb)smb_core(femmodel);

		if(ismasstransport){
			masstransport_core(femmodel);
			femmodel->UpdateVertexPositionsx();
		}
		
		if(isgroundingline){
			if(VerboseSolution()) _printf0_("   computing new grounding line position\n");
			GroundinglineMigrationx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters);

			femmodel->parameters->SetParam(MaskGroundediceLevelsetEnum,InputToExtrudeEnum);
			extrudefrombase_core(femmodel);
			femmodel->parameters->SetParam(BaseEnum,InputToExtrudeEnum);
			extrudefrombase_core(femmodel);
			femmodel->parameters->SetParam(SurfaceEnum,InputToExtrudeEnum);
			extrudefrombase_core(femmodel);
				
			if(save_results){
				int outputs[3] = {SurfaceEnum,BaseEnum,MaskGroundediceLevelsetEnum};
				femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],3);
			}
		}

		/*Calculate new basal melting on floating ice*/
		FloatingiceMeltingRatex(femmodel);
		
		if(isgia){
			if(VerboseSolution()) _printf0_("   computing glacial isostatic adjustment\n");
			#ifdef _HAVE_GIAIVINS_
			gia_core(femmodel);
			#else
			_error_("ISSM was not compiled with gia capabilities. Exiting");
			#endif
		}

		/*Sea level rise: */
		if(isslr || iscoupler) sealevelrise_core(femmodel);

		/*unload results*/
		if(VerboseSolution()) _printf0_("   computing requested outputs\n");
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs,save_results);
		if(isgroundingline && (groundingline_migration==SubelementMigrationEnum || groundingline_migration==SubelementMigration2Enum)){
			int outputs[1] = {MaskGroundediceLevelsetEnum};
			femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],1,save_results);
		}
		
		if(save_results){
			if(VerboseSolution()) _printf0_("   saving temporary results\n");
			OutputResultsx(femmodel);
		}

		if(recording_frequency && (step%recording_frequency==0)){
			if(VerboseSolution()) _printf0_("   checkpointing model \n");
			femmodel->CheckPoint();
		}

		/*Adaptive mesh refinement*/
		#ifdef _HAVE_NEOPZ_
		if(amr_frequency){
			if(save_results) femmodel->WriteMeshInResults();
			if(step%amr_frequency==0 && time<finaltime) femmodel->ReMesh();//Do not refine the last step
		}
		#endif
	
	}
	
	femmodel->RequestedDependentsx();

	/*Free ressources:*/	
	if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}
}
