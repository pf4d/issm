/*!\file: CreateParameters.cpp
 * \brief general driver for creating parameters dataset
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "../MeshPartitionx/MeshPartitionx.h"
#include "../ParseToolkitsOptionsx/ParseToolkitsOptionsx.h"
#include "./ModelProcessorx.h"

void CreateParameters(Parameters* parameters,IoModel* iomodel,char* rootpath,FILE* toolkitsoptionsfid,const int solution_type){

	int         i,j,m,k;
	int         numoutputs,materialtype,smb_model,basalforcing_model;
	char**      requestedoutputs = NULL;
	IssmDouble  time;

	/*parameters for mass flux:*/
	int          mass_flux_num_profiles     = 0;
	bool         qmu_mass_flux_present      = false;
	bool         autodiff_mass_flux_present = false;
	bool         mass_flux_present          = false;
	IssmDouble **array                      = NULL;
	int         *mdims_array                = NULL;
	int         *ndims_array                = NULL;
	IssmDouble  *temp_matrix                = NULL;
	int          temp_m,temp_n;
	IssmDouble  *matrix                     = NULL;
	int          count;

	IssmDouble *temp = NULL;
	IssmDouble  yts;
	int         N,M;

	/*Copy some constants from iomodel */
	parameters->AddObject(iomodel->CopyConstantObject("md.mesh.domain_type",DomainTypeEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.mesh.domain_dimension",DomainDimensionEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.output_frequency",SettingsOutputFrequencyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.recording_frequency",SettingsRecordingFrequencyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.constants.yts",ConstantsYtsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.start_time",TimesteppingStartTimeEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.final_time",TimesteppingFinalTimeEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.time_adapt",TimesteppingTimeAdaptEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.time_step",TimesteppingTimeStepEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.cfl_coefficient",TimesteppingCflCoefficientEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.interp_forcings",TimesteppingInterpForcingsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.lowmem",SettingsLowmemEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.debug.profiling",DebugProfilingEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.mesh.average_vertex_connectivity",MeshAverageVertexConnectivityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.waitonlock",SettingsWaitonlockEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.mesh.numberofelements",MeshNumberofelementsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.mesh.numberofvertices",MeshNumberofverticesEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.results_on_nodes",SettingsResultsOnNodesEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.io_gather",SettingsIoGatherEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.solver_residue_threshold",SettingsSolverResidueThresholdEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.isautodiff",AutodiffIsautodiffEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.qmu.isdakota",QmuIsdakotaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.inversion.iscontrol",InversionIscontrolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.inversion.type",InversionTypeEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.calving.law",CalvingLawEnum));
	{/*This is specific to ice...*/
		parameters->AddObject(iomodel->CopyConstantObject("md.mesh.elementtype",MeshElementtypeEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.steadystate.reltol",SteadystateReltolEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.steadystate.maxiter",SteadystateMaxiterEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.constants.referencetemperature",ConstantsReferencetemperatureEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.groundingline.migration",GroundinglineMigrationEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isstressbalance",TransientIsstressbalanceEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.ismasstransport",TransientIsmasstransportEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.issmb",TransientIssmbEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isthermal",TransientIsthermalEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isgroundingline",TransientIsgroundinglineEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isgia",TransientIsgiaEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isesa",TransientIsesaEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isdamageevolution",TransientIsdamageevolutionEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.ishydrology",TransientIshydrologyEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.ismovingfront",TransientIsmovingfrontEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isslr",TransientIsslrEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.iscoupler",TransientIscouplerEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.amr_frequency",TransientAmrFrequencyEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.materials.rheology_law",MaterialsRheologyLawEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.gia.cross_section_shape",GiaCrossSectionShapeEnum));
		/*amr properties*/	
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.level_max",AmrLevelMaxEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.region_level_1",AmrRegionLevel1Enum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.region_level_max",AmrRegionLevelMaxEnum));

		/*For stress balance only*/
		parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isFS",FlowequationIsFSEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.rift_penalty_threshold",StressbalanceRiftPenaltyThresholdEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.maxiter",StressbalanceMaxiterEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.restol",StressbalanceRestolEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.reltol",StressbalanceReltolEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.abstol",StressbalanceAbstolEnum));
		if(iomodel->domaintype==Domain3DEnum)
		 parameters->AddObject(iomodel->CopyConstantObject("md.mesh.numberoflayers",MeshNumberoflayersEnum));
	}

	
	/*Basal forcing parameters*/
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.model",BasalforcingsEnum));
	iomodel->FindConstant(&basalforcing_model,"md.basalforcings.model");
	switch(basalforcing_model){
		case FloatingMeltRateEnum:
			/*Nothing to add to parameters*/
			break;
		case LinearFloatingMeltRateEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.deepwater_melting_rate",BasalforcingsDeepwaterMeltingRateEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.deepwater_elevation",BasalforcingsDeepwaterElevationEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.upperwater_elevation",BasalforcingsUpperwaterElevationEnum));
			break;
		case MismipFloatingMeltRateEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.meltrate_factor",BasalforcingsMeltrateFactorEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.threshold_thickness",BasalforcingsThresholdThicknessEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.upperdepth_melt",BasalforcingsUpperdepthMeltEnum));
			break;
		case MantlePlumeGeothermalFluxEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.mantleconductivity",BasalforcingsMantleconductivityEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.nusselt",BasalforcingsNusseltEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.dtbg",BasalforcingsDtbgEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.plumeradius",BasalforcingsPlumeradiusEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.topplumedepth",BasalforcingsTopplumedepthEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.bottomplumedepth",BasalforcingsBottomplumedepthEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.plumex",BasalforcingsPlumexEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.plumey",BasalforcingsPlumeyEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.crustthickness",BasalforcingsCrustthicknessEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.uppercrustthickness",BasalforcingsUppercrustthicknessEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.uppercrustheat",BasalforcingsUppercrustheatEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.lowercrustheat",BasalforcingsLowercrustheatEnum));
			break;
		default:
			_error_("Basal forcing model "<<EnumToStringx(smb_model)<<" not supported yet");
	}

	/*some parameters that did not come with the iomodel: */
	parameters->AddObject(new IntParam(SolutionTypeEnum,solution_type));

	iomodel->FindConstant(&time,"md.timestepping.start_time");
	parameters->AddObject(new DoubleParam(TimeEnum,time));  
	parameters->AddObject(new IntParam(StepEnum,0));  

	/*By default, save all results*/
	parameters->AddObject(new BoolParam(SaveResultsEnum,true));

	/*Requested outputs */
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.transient.requested_outputs");
	parameters->AddObject(new IntParam(TransientNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(TransientRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.transient.requested_outputs");

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.steadystate.requested_outputs");
	parameters->AddObject(new IntParam(SteadystateNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(SteadystateRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.steadystate.requested_outputs");

	iomodel->FindConstant(&materialtype,"md.materials.type");
	if(materialtype==MatdamageiceEnum){
		iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.damage.requested_outputs");
		parameters->AddObject(new IntParam(DamageEvolutionNumRequestedOutputsEnum,numoutputs));
		if(numoutputs)parameters->AddObject(new StringArrayParam(DamageEvolutionRequestedOutputsEnum,requestedoutputs,numoutputs));
		iomodel->DeleteData(&requestedoutputs,numoutputs,"md.damage.requested_outputs");
	}

	/*Deal with mass flux segments: {{{*/
	iomodel->FetchData(&qmu_mass_flux_present,"md.qmu.mass_flux_segments_present");
	iomodel->FetchData(&autodiff_mass_flux_present,"md.autodiff.mass_flux_segments_present");

	if(qmu_mass_flux_present || autodiff_mass_flux_present)mass_flux_present=true;
	else mass_flux_present=false;
	parameters->AddObject(new BoolParam(MassFluxSegmentsPresentEnum,mass_flux_present));

	if(mass_flux_present){

		/*Fetch the mass flux segments necessary to compute the mass fluxes.  Build a DoubleMatArrayParam object out of them: */ 
		iomodel->FetchData(&array,&mdims_array,&ndims_array,&mass_flux_num_profiles,"md.qmu.mass_flux_segments");
		if(mass_flux_num_profiles==0)_error_("mass_flux_num_profiles is 0, when MassFlux computations were requested!");

		/*Go through segments, and extract those that belong to this cpu: */
		for(i=0;i<mass_flux_num_profiles;i++){
			temp_matrix=array[i];
			temp_m=mdims_array[i];
			temp_n=ndims_array[i];
			_assert_(temp_n==5);

			m=0;
			for(j=0;j<temp_m;j++){
				if (  iomodel->my_elements[reCast<int>(*(temp_matrix+5*j+4))-1] )m++;
			}
			if(m){
				matrix=xNewZeroInit<IssmDouble>(5*m);
				count=0;
				for(j=0;j<temp_m;j++){
					if (iomodel->my_elements[reCast<int>(*(temp_matrix+5*j+4))-1]){
						for(k=0;k<5;k++)*(matrix+5*count+k)=*(temp_matrix+5*j+k);
						count++;
					}
				}
			}
			else{
				matrix=NULL;
			}

			/*Assign: */
			array[i]=matrix;
			mdims_array[i]=m;
			ndims_array[i]=5;

			/*Free temporary matrix: */
			xDelete<IssmDouble>(temp_matrix);
		}

		/*Ok, we have an array of segments, different on every cpu. Create a DoubleMatArrayParam object with it: */
		parameters->AddObject(new DoubleMatArrayParam(MassFluxSegmentsEnum,array,mass_flux_num_profiles,mdims_array,ndims_array));

		/*Free data: */
		for(i=0;i<mass_flux_num_profiles;i++){
			IssmDouble* matrix=array[i];
			xDelete<IssmDouble>(matrix);
		}
		xDelete<int>(mdims_array); 
		xDelete<int>(ndims_array);
		xDelete<IssmDouble*>(array);
	}
	/*}}}*/

	/*Before returning, create parameters in case we are running Qmu or control types runs: */
	CreateParametersControl(parameters,iomodel,solution_type);

	#ifdef _HAVE_DAKOTA_
	CreateParametersDakota(parameters,iomodel,rootpath);
	#endif

	/*Now, deal with toolkits options, which need to be put into the parameters dataset: */
	ParseToolkitsOptionsx(parameters,toolkitsoptionsfid);

 	#ifdef _HAVE_ADOLC_
	CreateParametersAutodiff(parameters,iomodel);
	#endif
}
