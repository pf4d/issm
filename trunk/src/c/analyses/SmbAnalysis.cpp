#include "./SmbAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void SmbAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*No constraints*/
}/*}}}*/
void SmbAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*No loads*/
}/*}}}*/
void SmbAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/
	::CreateNodes(nodes,iomodel,SmbAnalysisEnum,P1Enum);
}/*}}}*/
int  SmbAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void SmbAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/
	
	int    smb_model;
	bool   isdelta18o,ismungsm,isd18opd;
	
	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}
	
	/*Figure out smb model: */
	iomodel->FindConstant(&smb_model,"md.smb.model");
			
	switch(smb_model){
		case SMBforcingEnum:
			iomodel->FetchDataToInput(elements,"md.smb.mass_balance",SmbMassBalanceEnum,0.);
			break;
		case SMBgembEnum:
			iomodel->FetchDataToInput(elements,"md.smb.Ta",SmbTaEnum);
			iomodel->FetchDataToInput(elements,"md.smb.V",SmbVEnum);
			iomodel->FetchDataToInput(elements,"md.smb.dswrf",SmbDswrfEnum);
			iomodel->FetchDataToInput(elements,"md.smb.dlwrf",SmbDlwrfEnum);
			iomodel->FetchDataToInput(elements,"md.smb.P",SmbPEnum);
			iomodel->FetchDataToInput(elements,"md.smb.eAir",SmbEAirEnum);
			iomodel->FetchDataToInput(elements,"md.smb.pAir",SmbPAirEnum);
			iomodel->FetchDataToInput(elements,"md.smb.zTop",SmbZTopEnum);
			iomodel->FetchDataToInput(elements,"md.smb.dzTop",SmbDzTopEnum);
			iomodel->FetchDataToInput(elements,"md.smb.dzMin",SmbDzMinEnum);
			iomodel->FetchDataToInput(elements,"md.smb.zY",SmbZYEnum);
			iomodel->FetchDataToInput(elements,"md.smb.zMax",SmbZMaxEnum);
			iomodel->FetchDataToInput(elements,"md.smb.zMin",SmbZMinEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Tmean",SmbTmeanEnum);
			iomodel->FetchDataToInput(elements,"md.smb.C",SmbCEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Tz",SmbTzEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Vz",SmbVzEnum);
			InputUpdateFromConstantx(elements,0.,SmbIsInitializedEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Dzini",SmbDziniEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Dini",SmbDiniEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Reini",SmbReiniEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Gdnini",SmbGdniniEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Gspini",SmbGspiniEnum);
			iomodel->FetchDataToInput(elements,"md.smb.ECini",SmbECiniEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Wini",SmbWiniEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Aini",SmbAiniEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Tini",SmbTiniEnum);
			iomodel->FetchDataToInput(elements,"md.smb.Sizeini",SmbSizeiniEnum);
			break;
		case SMBpddEnum:
			iomodel->FindConstant(&isdelta18o,"md.smb.isdelta18o");
			iomodel->FindConstant(&ismungsm,"md.smb.ismungsm");
			iomodel->FetchDataToInput(elements,"md.thermal.spctemperature",ThermalSpctemperatureEnum);
			iomodel->FetchDataToInput(elements,"md.smb.s0p",SmbS0pEnum);
			iomodel->FetchDataToInput(elements,"md.smb.s0t",SmbS0tEnum);
			if(isdelta18o || ismungsm){
				iomodel->FetchDataToInput(elements,"md.smb.temperatures_lgm",SmbTemperaturesLgmEnum);
				iomodel->FetchDataToInput(elements,"md.smb.temperatures_presentday",SmbTemperaturesPresentdayEnum);
				iomodel->FetchDataToInput(elements,"md.smb.precipitations_presentday",SmbPrecipitationsPresentdayEnum);
				iomodel->FetchDataToInput(elements,"md.smb.precipitations_lgm",SmbPrecipitationsLgmEnum);
			}
			else{
				iomodel->FetchDataToInput(elements,"md.smb.precipitation",SmbPrecipitationEnum);
				iomodel->FetchDataToInput(elements,"md.smb.monthlytemperatures",SmbMonthlytemperaturesEnum);
			}
			break;
		case SMBd18opddEnum:
			iomodel->FindConstant(&ismungsm,"md.smb.ismungsm");
			iomodel->FindConstant(&isd18opd,"md.smb.isd18opd");
			iomodel->FetchDataToInput(elements,"md.thermal.spctemperature",ThermalSpctemperatureEnum);
			iomodel->FetchDataToInput(elements,"md.smb.s0p",SmbS0pEnum);
			iomodel->FetchDataToInput(elements,"md.smb.s0t",SmbS0tEnum);
			if(isd18opd){
				iomodel->FetchDataToInput(elements,"md.smb.temperatures_presentday",SmbTemperaturesPresentdayEnum);
				iomodel->FetchDataToInput(elements,"md.smb.precipitations_presentday",SmbPrecipitationsPresentdayEnum);
			}

			break;
		case SMBgradientsEnum:
			iomodel->FetchDataToInput(elements,"md.smb.href",SmbHrefEnum);
			iomodel->FetchDataToInput(elements,"md.smb.smbref",SmbSmbrefEnum);
			iomodel->FetchDataToInput(elements,"md.smb.b_pos",SmbBPosEnum);
			iomodel->FetchDataToInput(elements,"md.smb.b_neg",SmbBNegEnum);
			break;
		case SMBgradientselaEnum:
			iomodel->FetchDataToInput(elements,"md.smb.ela",SmbElaEnum);
			iomodel->FetchDataToInput(elements,"md.smb.b_pos",SmbBPosEnum);
			iomodel->FetchDataToInput(elements,"md.smb.b_neg",SmbBNegEnum);
			iomodel->FetchDataToInput(elements,"md.smb.b_max",SmbBMaxEnum);
			iomodel->FetchDataToInput(elements,"md.smb.b_min",SmbBMinEnum);
			break;
		case SMBhenningEnum:
			iomodel->FetchDataToInput(elements,"md.smb.smbref",SmbSmbrefEnum,0.);
			break;
		case SMBcomponentsEnum:
			iomodel->FetchDataToInput(elements,"md.smb.accumulation",SmbAccumulationEnum,0.);
			iomodel->FetchDataToInput(elements,"md.smb.evaporation",SmbEvaporationEnum,0.);
			iomodel->FetchDataToInput(elements,"md.smb.runoff",SmbRunoffEnum,0.);
			break;
		case SMBmeltcomponentsEnum:
			iomodel->FetchDataToInput(elements,"md.smb.accumulation",SmbAccumulationEnum,0.);
			iomodel->FetchDataToInput(elements,"md.smb.evaporation",SmbEvaporationEnum,0.);
			iomodel->FetchDataToInput(elements,"md.smb.melt",SmbMeltEnum,0.);
			iomodel->FetchDataToInput(elements,"md.smb.refreeze",SmbRefreezeEnum,0.);
			break;
		default:
			_error_("Surface mass balance model "<<EnumToStringx(smb_model)<<" not supported yet");
	}
	
	

}/*}}}*/
void SmbAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;
	bool    isdelta18o,ismungsm,isd18opd,interp;
	int     smb_model;
	IssmDouble *temp = NULL;
	int         N,M;
	
	parameters->AddObject(iomodel->CopyConstantObject("md.smb.model",SmbEnum));
	
	iomodel->FindConstant(&smb_model,"md.smb.model");
	iomodel->FindConstant(&interp,"md.timestepping.interp_forcings");
	
	switch(smb_model){
		case SMBforcingEnum:
			/*Nothing to add to parameters*/
			break;
		case SMBgembEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.aIdx",SmbAIdxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.swIdx",SmbSwIdxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.denIdx",SmbDenIdxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.outputFreq",SmbOutputFreqEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.cldFrac",SmbCldFracEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.t0wet",SmbT0wetEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.t0dry",SmbT0dryEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.K",SmbKEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.aSnow",SmbASnowEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.aIce",SmbAIceEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.dt",SmbDtEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isgraingrowth",SmbIsgraingrowthEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isalbedo",SmbIsalbedoEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isshortwave",SmbIsshortwaveEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isthermal",SmbIsthermalEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isaccumulation",SmbIsaccumulationEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ismelt",SmbIsmeltEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isdensification",SmbIsdensificationEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isturbulentflux",SmbIsturbulentfluxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.InitDensityScaling",SmbInitDensityScalingEnum));
			break;
		case SMBpddEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isdelta18o",SmbIsdelta18oEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ismungsm",SmbIsmungsmEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.desfac",SmbDesfacEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlaps",SmbRlapsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlapslgm",SmbRlapslgmEnum));
			iomodel->FindConstant(&isdelta18o,"md.smb.isdelta18o");
			iomodel->FindConstant(&ismungsm,"md.smb.ismungsm");

			if(ismungsm){
			  iomodel->FetchData(&temp,&N,&M,"md.smb.Pfac"); _assert_(N==2);
			  parameters->AddObject(new TransientParam(SmbPfacEnum,&temp[0],&temp[M],interp,M));
			  iomodel->DeleteData(temp,"md.smb.Pfac");
			
			  iomodel->FetchData(&temp,&N,&M,"md.smb.Tdiff"); _assert_(N==2);
			  parameters->AddObject(new TransientParam(SmbTdiffEnum,&temp[0],&temp[M],interp,M));
			  iomodel->DeleteData(temp,"md.smb.Tdiff");

			  iomodel->FetchData(&temp,&N,&M,"md.smb.sealev"); _assert_(N==2);
			  parameters->AddObject(new TransientParam(SmbSealevEnum,&temp[0],&temp[M],interp,M));
			  iomodel->DeleteData(temp,"md.smb.sealev");
			}
			if(isdelta18o){
				iomodel->FetchData(&temp,&N,&M,"md.smb.delta18o"); _assert_(N==2);
				parameters->AddObject(new TransientParam(SmbDelta18oEnum,&temp[0],&temp[M],interp,M));
				iomodel->DeleteData(temp,"md.smb.delta18o");

				iomodel->FetchData(&temp,&N,&M,"md.smb.delta18o_surface"); _assert_(N==2);
				parameters->AddObject(new TransientParam(SmbDelta18oSurfaceEnum,&temp[0],&temp[M],interp,M));
				iomodel->DeleteData(temp,"md.smb.delta18o_surface");
			}
			break;
		case SMBd18opddEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ismungsm",SmbIsmungsmEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isd18opd",SmbIsd18opdEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.desfac",SmbDesfacEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlaps",SmbRlapsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlapslgm",SmbRlapslgmEnum));
			iomodel->FindConstant(&ismungsm,"md.smb.ismungsm");
			iomodel->FindConstant(&isd18opd,"md.smb.isd18opd");
			if(isd18opd){
				iomodel->FetchData(&temp,&N,&M,"md.smb.delta18o"); _assert_(N==2);
				parameters->AddObject(new TransientParam(SmbDelta18oEnum,&temp[0],&temp[M],interp,M));
				iomodel->DeleteData(temp,"md.smb.delta18o");
				
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.dpermil",SmbDpermilEnum));
			   parameters->AddObject(iomodel->CopyConstantObject("md.smb.f",SmbFEnum));
			}
			break;
		case SMBgradientsEnum:
			/*Nothing to add to parameters*/
			break;
		case SMBgradientselaEnum:
			/*Nothing to add to parameters*/
			break;
		case SMBhenningEnum:
			/*Nothing to add to parameters*/
			break;
		case SMBcomponentsEnum:
			/*Nothing to add to parameters*/
			break;
		case SMBmeltcomponentsEnum:
			/*Nothing to add to parameters*/
			break;
		default:
			_error_("Surface mass balance model "<<EnumToStringx(smb_model)<<" not supported yet");
	}

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.smb.requested_outputs");
	parameters->AddObject(new IntParam(SmbNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(SmbRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.smb.requested_outputs");

}/*}}}*/

/*Finite Element Analysis*/
void           SmbAnalysis::Core(FemModel* femmodel){/*{{{*/

	int    smb_model;

	/*Figure out smb model: */
	femmodel->parameters->FindParam(&smb_model,SmbEnum);
	
	/*branch to correct module*/
	switch(smb_model){
		case SMBforcingEnum:
			/*Nothing to be done*/
			break;
		case SMBgembEnum:
			Gembx(femmodel);
			break;
		case SMBpddEnum:
			bool isdelta18o,ismungsm;
			femmodel->parameters->FindParam(&isdelta18o,SmbIsdelta18oEnum);
			femmodel->parameters->FindParam(&ismungsm,SmbIsmungsmEnum);
			if(isdelta18o){
				if(VerboseSolution()) _printf0_("   call Delta18oParameterization module\n");
				Delta18oParameterizationx(femmodel);
			} 
			if(ismungsm){
				if(VerboseSolution()) _printf0_("   call MungsmtpParameterization module\n");
				MungsmtpParameterizationx(femmodel);
			} 
			if(VerboseSolution()) _printf0_("   call positive degree day module\n");
			PositiveDegreeDayx(femmodel);
			break;
		case SMBd18opddEnum:
			bool isd18opd;
			femmodel->parameters->FindParam(&isd18opd,SmbIsd18opdEnum);
			if(isd18opd){
				if(VerboseSolution()) _printf0_("   call Delta18opdParameterization module\n");
				Delta18opdParameterizationx(femmodel);
				if(VerboseSolution()) _printf0_("   call positive degree day module\n");
				PositiveDegreeDayx(femmodel);
			} 
			break;
		case SMBgradientsEnum:
			if(VerboseSolution())_printf0_("	call smb gradients module\n");
			SmbGradientsx(femmodel);
			break;
		case SMBgradientselaEnum:
			if(VerboseSolution())_printf0_("	call smb gradients ela module\n");
			SmbGradientsElax(femmodel);
			break;
		case SMBhenningEnum:
			if(VerboseSolution())_printf0_("  call smb Henning module\n");
			SmbHenningx(femmodel);
			break;
		case SMBcomponentsEnum:
			if(VerboseSolution())_printf0_("  call smb Components module\n");
			SmbComponentsx(femmodel);
			break;
		case SMBmeltcomponentsEnum:
			if(VerboseSolution())_printf0_("  call smb Melt Components module\n");
			SmbMeltComponentsx(femmodel);
			break;
		case SMBgcmEnum:
			/*Nothing to be done*/
			break;
		default:
			_error_("Surface mass balance model "<<EnumToStringx(smb_model)<<" not supported yet");
	}

}/*}}}*/
ElementVector* SmbAnalysis::CreateDVector(Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementMatrix* SmbAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* SmbAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* SmbAnalysis::CreatePVector(Element* element){/*{{{*/
_error_("not implemented yet");
}/*}}}*/
void           SmbAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           SmbAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           SmbAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           SmbAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
