/*!\file: levelset_core.cpp
 * \brief: levelset-module to update the ice domain
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../solutionsequences/solutionsequences.h"
#include "../modules/modules.h"

void movingfront_core(FemModel* femmodel){

	/* intermediaries */
	bool save_results,isstressbalance,ismasstransport,isthermal,isenthalpy,islevelset,ismovingfront;
	int  domaintype, num_extrapol_vars, index,reinit_frequency,step;
	int* extrapol_vars=NULL;
	Analysis  *analysis=NULL;

	/* recover parameters */
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&isstressbalance,TransientIsstressbalanceEnum);
	femmodel->parameters->FindParam(&ismasstransport,TransientIsmasstransportEnum);
	femmodel->parameters->FindParam(&isthermal,TransientIsthermalEnum);
	femmodel->parameters->FindParam(&ismovingfront,TransientIsmovingfrontEnum);
	femmodel->parameters->FindParam(&reinit_frequency,LevelsetReinitFrequencyEnum);
	femmodel->parameters->FindParam(&step,StepEnum);
	if(isthermal && domaintype==Domain3DEnum) femmodel->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);

	if(!ismovingfront) return;

	/* start the work from here */
	Calvingx(femmodel);
	if(VerboseSolution()) _printf0_("   computing level set transport\n");

	/* smoothen slope of lsf for computation of normal on ice domain*/
	levelsetfunctionslope_core(femmodel);

	/* determine variables for extrapolation */
	num_extrapol_vars=0;
	if(isstressbalance) num_extrapol_vars+=2;
	if(ismasstransport) num_extrapol_vars+=1;
	if(isthermal && domaintype==Domain3DEnum) num_extrapol_vars+=1;
	extrapol_vars=xNew<int>(num_extrapol_vars);
	index=0;
	if(isstressbalance){
		extrapol_vars[index]=VxEnum; index++;
		extrapol_vars[index]=VyEnum; index++;
	}
	if(ismasstransport){extrapol_vars[index]=ThicknessEnum; index++;}
	if(isthermal && domaintype==Domain3DEnum){
		if(isenthalpy){extrapol_vars[index]=EnthalpyEnum;}
		else{extrapol_vars[index]=TemperatureEnum;}
		index++;
	}

	/* extrapolate */
	analysis = new ExtrapolationAnalysis();
	for(int iv=0;iv<num_extrapol_vars;iv++){
		femmodel->parameters->SetParam(extrapol_vars[iv],ExtrapolationVariableEnum); 
		analysis->Core(femmodel);
	}
	xDelete<int>(extrapol_vars);
	delete analysis;	

	/* solve level set equation */
	analysis = new LevelsetAnalysis();
	analysis->Core(femmodel);
	delete analysis;

	/*Reset levelset if needed*/
	if(reinit_frequency && (step%reinit_frequency==0)){
		if(VerboseSolution()) _printf0_("   reinitializing level set\n");
		femmodel->ResetLevelset();
	}

	/* update vertices included for next calculation */
	GetMaskOfIceVerticesLSMx(femmodel);

	/* add computation domain mask to outputs */
	if(save_results){
		int outputs[1] = {IceMaskNodeActivationEnum};
		femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],1);
	}
}
