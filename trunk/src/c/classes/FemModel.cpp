/*!\file FemModel.cpp
 * \brief: implementation of the FemModel object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include "../cores/cores.h"
#include "../shared/io/io.h"
#include "./classes.h"
#include "./modules/modules.h"
#include "../shared/Enum/Enum.h"
#include "../analyses/analyses.h"

/*module includes: {{{*/
#include "../modules/ModelProcessorx/ModelProcessorx.h"
#include "../modules/VerticesDofx/VerticesDofx.h"
#include "../modules/SpcNodesx/SpcNodesx.h"
#include "../modules/ConfigureObjectsx/ConfigureObjectsx.h"
#include "../modules/ParseToolkitsOptionsx/ParseToolkitsOptionsx.h"
#include "../modules/GetVectorFromInputsx/GetVectorFromInputsx.h"
#include "../modules/InputUpdateFromVectorx/InputUpdateFromVectorx.h"
#include "../modules/NodesDofx/NodesDofx.h"
#include "../modules/SurfaceAbsVelMisfitx/SurfaceAbsVelMisfitx.h"
#include "../modules/SurfaceRelVelMisfitx/SurfaceRelVelMisfitx.h"
#include "../modules/SurfaceLogVelMisfitx/SurfaceLogVelMisfitx.h"
#include "../modules/SurfaceLogVxVyMisfitx/SurfaceLogVxVyMisfitx.h"
#include "../modules/SurfaceAverageVelMisfitx/SurfaceAverageVelMisfitx.h"
#include "../modules/ThicknessAbsMisfitx/ThicknessAbsMisfitx.h"
#include "../modules/ThicknessAlongGradientx/ThicknessAlongGradientx.h"
#include "../modules/ThicknessAcrossGradientx/ThicknessAcrossGradientx.h"
#include "../modules/RheologyBbarAbsGradientx/RheologyBbarAbsGradientx.h"
#include "../modules/DragCoefficientAbsGradientx/DragCoefficientAbsGradientx.h"
#include "../modules/NodalValuex/NodalValuex.h"
#include "../modules/GetVectorFromInputsx/GetVectorFromInputsx.h"
#include "../modules/AverageOntoPartitionx/AverageOntoPartitionx.h"
/*}}}*/

/*Object constructors and destructor*/
FemModel::FemModel(int argc,char** argv,ISSM_MPI_Comm incomm,bool trace){/*{{{*/

	/*configuration: */
	int  solution_type;
	int  ierr;

	/*File names*/
	char *lockfilename   = NULL;
	char *binfilename    = NULL;
	char *outbinfilename = NULL;
	char *petscfilename  = NULL;
	char *restartfilename  = NULL;
	char *rootpath       = NULL;

	/*First things first, store the communicator, and set it as a global variable: */
	IssmComm::SetComm(incomm);

	/*Now, initialize PETSC: */
	#ifdef _HAVE_PETSC_
	PETSC_COMM_WORLD=incomm;
	ierr=PetscInitialize(&argc,&argv,(char*)0,"");  if(ierr) _error_("Could not initialize Petsc");
	#endif

	/*Start profiler: */
	this->profiler=new Profiler();
	profiler->Tag(START);

	/*From command line arguments, retrieve different filenames needed to create the FemModel: */
	ProcessArguments(&solution_type,&binfilename,&outbinfilename,&petscfilename,&lockfilename,&restartfilename,&rootpath,argc,argv);

	/*Create femmodel from input files: */
	profiler->Tag(STARTINIT);
	this->InitFromFiles(rootpath,binfilename,outbinfilename,petscfilename,lockfilename,restartfilename, solution_type,trace,NULL);
	profiler->Tag(FINISHINIT);

	/*Save communicator in the parameters dataset: */
	this->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(incomm,FemModelCommEnum));

	#ifdef _HAVE_NEOPZ_
	this->InitializeAdaptiveRefinement();
	#endif

	/*Free resources */
	xDelete<char>(lockfilename);
	xDelete<char>(binfilename);
	xDelete<char>(outbinfilename);
	xDelete<char>(petscfilename);
	xDelete<char>(restartfilename);
	xDelete<char>(rootpath);

}
/*}}}*/
FemModel::FemModel(char* rootpath, char* inputfilename, char* outputfilename, char* toolkitsfilename, char* lockfilename, char* restartfilename, ISSM_MPI_Comm incomm, int solution_type,IssmPDouble* X){ /*{{{*/

	bool traceon=true;
	this->profiler=NULL; /*avoid leak, as we are not using the profiler ever in ad control run. */
	
	/*Store the communicator, but do not set it as a global variable, as this has already 
	 * been done by the FemModel that called this copy constructor: */
	IssmComm::SetComm(incomm);

	/*Create femmodel from input files, with trace activated: */
	profiler->Tag(STARTINIT);
	this->InitFromFiles(rootpath,inputfilename,outputfilename,toolkitsfilename,lockfilename,restartfilename, solution_type,traceon,X);
	profiler->Tag(FINISHINIT);
	
	/*Save communicator in the parameters dataset: */
	this->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(incomm,FemModelCommEnum));

}
/*}}}*/
FemModel::~FemModel(){/*{{{*/

	/*Intermediary*/
	FILE *output_fid;
	char *outbinfilename = NULL;
	char *lockfilename   = NULL;

	#ifndef _HAVE_JAVASCRIPT_
	if(this->parameters->Exist(OutputFileNameEnum)) this->parameters->FindParam(&outbinfilename,OutputFileNameEnum);
	if(this-parameters->Exist(LockFileNameEnum)) this->parameters->FindParam(&lockfilename,LockFileNameEnum);
	#endif

	/*Delete all the datasets: */
	if(analysis_type_list)xDelete<int>(analysis_type_list);
	if(outbinfilename)xDelete<char>(outbinfilename);
	if(lockfilename)xDelete<char>(lockfilename);
	if(elements)delete elements;
	if(nodes)delete nodes;
	if(vertices)delete vertices;
	if(constraints)delete constraints;
	if(loads)delete loads;
	if(materials)delete materials;
	if(parameters)delete parameters;
	if(results)delete results;
	
	#ifdef _HAVE_NEOPZ_
	if(amr)delete amr;
	#endif

	/*Now delete: */
	if(profiler)delete profiler;
	
	
}
/*}}}*/

/*Object management*/
void FemModel::CheckPoint(void){/*{{{*/

	FILE* restartfid=NULL;
	char* restartfilename = NULL;
	int   femmodel_size;
	char* femmodel_buffer=NULL;
	char* femmodel_buffer_ini=NULL;

	/*First, recover the name of the restart file: */
	parameters->FindParam(&restartfilename,RestartFileNameEnum);
	
	/*Open file for writing: */
	restartfid=pfopen(restartfilename,"wb");

	/*Initialize: */
	femmodel_size=0;

	/*Create buffer to hold marshalled femmodel: */
	this->Marshall(NULL,&femmodel_size,MARSHALLING_SIZE);
	femmodel_buffer=xNew<char>(femmodel_size); 
	/*Keep track of initial position of femmodel_buffer: */
	femmodel_buffer_ini=femmodel_buffer;
	
	/*Marshall:*/
	this->Marshall(&femmodel_buffer,NULL,MARSHALLING_FORWARD);

	/*Reset position of buffer: */
	femmodel_buffer=femmodel_buffer_ini;

	/*write buffer: */
	fwrite(femmodel_buffer,femmodel_size,sizeof(char),restartfid);

	/*Done, close file :*/
	pfclose(restartfid,restartfilename);

	/*Free ressources: */
	xDelete<char>(femmodel_buffer);
	xDelete<char>(restartfilename);

}
/*}}}*/
void FemModel::CleanUp(void){/*{{{*/

	/*Intermediary*/
	char *lockfilename   = NULL;
	bool  waitonlock     = false;


	/*Write lock file if requested: */
	this->parameters->FindParam(&waitonlock,SettingsWaitonlockEnum);
	this->parameters->FindParam(&lockfilename,LockFileNameEnum);
	if(waitonlock){
		_printf0_("write lock file:\n");
		WriteLockFile(lockfilename);
	}

	/*Before we delete the profiler, report statistics for this run: */
	profiler->Tag(FINISH);  //final tagging
	_printf0_("\n");
	_printf0_("   "<<setw(40)<<left<<"FemModel initialization elapsed time:"<<profiler->DeltaTime(STARTINIT,FINISHINIT) << "\n");
	_printf0_("   "<<setw(40)<<left<<"Core solution elapsed time:"<<profiler->DeltaTime(STARTCORE,FINISHCORE) << "\n");
	_printf0_("\n");
	_printf0_("   Total elapsed time: "
				<<profiler->DeltaTimeModHour(START,FINISH)<<" hrs "
				<<profiler->DeltaTimeModMin(START,FINISH)<<" min "
				<<profiler->DeltaTimeModSec(START,FINISH)<<" sec"
				);
	_printf0_("\n");

	/*Finalize PETSC for this model: */
	#ifdef _HAVE_PETSC_
	_printf0_("closing PETSc\n");
	PetscFinalize();
	#endif

	/*Clean up*/
	xDelete<char>(lockfilename);
} /*}}}*/
FemModel* FemModel::copy(void){/*{{{*/

	FemModel* output=NULL;
	int       i;
	int       analysis_type;

	output=new FemModel(*this); //Use default copy constructor.

	output->nummodels = this->nummodels;
	output->solution_type = this->solution_type;
	output->analysis_counter = this->analysis_counter;

	/*Now, deep copy arrays: */
	output->analysis_type_list=xNew<int>(nummodels);
	xMemCpy<int>(output->analysis_type_list,this->analysis_type_list,this->nummodels);

	output->profiler=static_cast<Profiler*>(this->profiler->copy());

	output->loads=static_cast<Loads*>(this->loads->Copy());
	output->materials=static_cast<Materials*>(this->materials->Copy());
	output->parameters=static_cast<Parameters*>(this->parameters->Copy());
	output->constraints=static_cast<Constraints*>(this->constraints->Copy());
	output->results=static_cast<Results*>(this->results->Copy());

	output->nodes=static_cast<Nodes*>(this->nodes->Copy());
	output->vertices=static_cast<Vertices*>(this->vertices->Copy());
	output->elements=static_cast<Elements*>(this->elements->Copy());

	/*reset hooks for elements, loads and nodes: */
	output->elements->ResetHooks();
	output->loads->ResetHooks();
	output->materials->ResetHooks();

	/*do the post-processing of the datasets to get an FemModel that can actually run analyses: */
	for(i=0;i<nummodels;i++){
		analysis_type=output->analysis_type_list[i];
		output->SetCurrentConfiguration(analysis_type);
		if(i==0) VerticesDofx(output->vertices,output->parameters); //only call once, we only have one set of vertices
		SpcNodesx(output->nodes,output->constraints,output->parameters,analysis_type);
		NodesDofx(output->nodes,output->parameters,analysis_type);
		ConfigureObjectsx(output->elements,output->loads,output->nodes,output->vertices,output->materials,output->parameters);
	}

	/*Reset current configuration: */
	analysis_type=output->analysis_type_list[analysis_counter];
	output->SetCurrentConfiguration(analysis_type);

	return output;
}
/*}}}*/
void FemModel::Echo(void){/*{{{*/

	_printf_("FemModel echo: \n");
	_printf_("   number of fem models: " << nummodels << "\n");
	_printf_("   analysis_type_list: \n");
	for(int i=0;i<nummodels;i++)_printf_("     " << i << ": " << EnumToStringx(analysis_type_list[i]) << "\n");
	_printf_("   current analysis_type: \n");
	_printf_("     " << analysis_counter << ": " << EnumToStringx(analysis_type_list[analysis_counter]) << "\n");

}
/*}}}*/
void FemModel::InitFromFiles(char* rootpath, char* inputfilename, char* outputfilename, char* toolkitsfilename, char* lockfilename, char* restartfilename, const int in_solution_type,bool trace,IssmPDouble* X){/*{{{*/

	/*intermediary*/
	FILE *IOMODEL            = NULL;
	FILE *toolkitsoptionsfid = NULL;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Open input file descriptor on cpu 0: */
	if(my_rank==0) IOMODEL = pfopen0(inputfilename ,"rb");

	/*Open toolkits file: */
	toolkitsoptionsfid=pfopen(toolkitsfilename,"r");

	/*Now, go create FemModel:*/
	this->InitFromFids(rootpath,IOMODEL,toolkitsoptionsfid,in_solution_type,trace,X);

	/*Close input file and toolkits file descriptors: */
	if(my_rank==0) pfclose(IOMODEL,inputfilename);
	pfclose(toolkitsoptionsfid,toolkitsfilename);

	/*Now save all of these file names into parameters, you never know when you might need them: */
	this->parameters->AddObject(new StringParam(ToolkitsFileNameEnum,toolkitsfilename));
	this->parameters->AddObject(new StringParam(RootPathEnum,rootpath));
	this->parameters->AddObject(new StringParam(InputFileNameEnum,inputfilename));
	this->parameters->AddObject(new StringParam(OutputFileNameEnum,outputfilename));
	this->parameters->AddObject(new StringParam(LockFileNameEnum,lockfilename));
	this->parameters->AddObject(new StringParam(RestartFileNameEnum,restartfilename));

}/*}}}*/
void FemModel::InitFromFids(char* rootpath, FILE* IOMODEL, FILE* toolkitsoptionsfid, int in_solution_type, bool trace, IssmPDouble* X){/*{{{*/
	
	/*Initialize internal data: */
	this->solution_type    = in_solution_type;
	this->analysis_counter = nummodels-1;   //point to last analysis_type carried out.
	this->results          = new Results(); //not initialized by CreateDataSets
	
	/*create IoModel */
	IoModel* iomodel = new IoModel(IOMODEL,in_solution_type,trace,X);

	/*Figure out what analyses are activated for this solution*/
	SolutionAnalysesList(&this->analysis_type_list,&this->nummodels,iomodel,this->solution_type);

	/*create datasets for all analyses*/
	ModelProcessorx(&this->elements,&this->nodes,&this->vertices,&this->materials,&this->constraints,&this->loads,&this->parameters,iomodel,toolkitsoptionsfid,rootpath,this->solution_type,this->nummodels,this->analysis_type_list);

	/*do the post-processing of the datasets to get an FemModel that can actually run analyses: */
	for(int i=0;i<nummodels;i++){

		if(VerboseMProcessor()) _printf0_("   Processing finite element model of analysis " << EnumToStringx(analysis_type_list[i]) << ":\n");
		this->SetCurrentConfiguration(analysis_type_list[i]);

		if(VerboseMProcessor()) _printf0_("      configuring element and loads\n");
		ConfigureObjectsx(elements, loads, nodes, vertices, materials,parameters);
		
		if(i==0){
			if(VerboseMProcessor()) _printf0_("      creating vertex PIDs\n");
			VerticesDofx(vertices,parameters); 

			if(VerboseMProcessor()) _printf0_("      detecting active vertices\n");
			GetMaskOfIceVerticesLSMx(this);
		}

		if(VerboseMProcessor()) _printf0_("      resolving node constraints\n");
		SpcNodesx(nodes,constraints,parameters,analysis_type_list[i]); 

		if(VerboseMProcessor()) _printf0_("      creating nodal degrees of freedom\n");
		NodesDofx(nodes,parameters,analysis_type_list[i]);
	}

	/*Clean up*/
	delete iomodel;
}/*}}}*/
void FemModel::Marshall(char** pmarshalled_data, int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	int       i;
	int       analysis_type;

	if(marshall_direction==MARSHALLING_BACKWARD){
		delete this->loads;
		delete this->materials;
		delete this->parameters;
		delete this->constraints;
		delete this->results;
		delete this->nodes;
		delete this->vertices;
		delete this->elements;
		xDelete<int>(this->analysis_type_list);

		this->loads       = new Loads();
		this->materials   = new Materials();
		this->parameters  = new Parameters();
		this->constraints = new Constraints();
		this->results     = new Results();
		this->nodes       = new Nodes();
		this->vertices    = new Vertices();
		this->elements    = new Elements();
	}

	MARSHALLING_ENUM(FemModelEnum);

	MARSHALLING(solution_type);
	MARSHALLING(analysis_counter);
	MARSHALLING(nummodels);
	MARSHALLING_DYNAMIC(analysis_type_list,int,nummodels);

	this->loads->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->materials->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->parameters->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->constraints->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->results->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->nodes->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->vertices->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	this->elements->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);

	if(marshall_direction==MARSHALLING_BACKWARD){
		/*reset hooks for elements, loads and nodes:*/
		this->elements->ResetHooks();
		this->loads->ResetHooks();
		this->materials->ResetHooks();

		/*do the post-processing of the datasets to get an FemModel that can actually run analyses:*/
		for(i=0;i<nummodels;i++){
			analysis_type=this->analysis_type_list[i];
			SetCurrentConfiguration(analysis_type);
			if(i==0) VerticesDofx(this->vertices,this->parameters); //only call once, we only have one set of vertices
			SpcNodesx(this->nodes,this->constraints,this->parameters,analysis_type);
			NodesDofx(this->nodes,this->parameters,analysis_type);
			ConfigureObjectsx(this->elements,this->loads,this->nodes,this->vertices,this->materials,this->parameters);
		}

		//Reset current configuration:
		analysis_type=this->analysis_type_list[analysis_counter];
		SetCurrentConfiguration(analysis_type);
	}

}
/*}}}*/
void FemModel::Restart(){ /*{{{*/

	FILE* restartfid=NULL;
	char* restartfilename = NULL;
	int   femmodel_size=0; 
	int   fread_return=0; 
	char* femmodel_buffer=NULL;
	char* femmodel_buffer_ini=NULL;

	/*First, recover the name of the restart file: */
	parameters->FindParam(&restartfilename,RestartFileNameEnum);

	/*Now, figure out whether this file actually exists!: */
	restartfid=pfopen(restartfilename,"r",false);

	if(restartfid==NULL){
		xDelete<char>(restartfilename);
		return; //could not find the file, so no restart possible.
	}

	/*Print banner*/
	_printf0_("                                                                    \n");
	_printf0_("====================================================================\n");
	_printf0_("                         RESTART DETECTED                           \n");
	_printf0_("                                                                    \n");
	_printf0_("  Restart file: "<<restartfilename<<"                               \n");
	_printf0_("====================================================================\n");
	_printf0_("                                                                    \n");

	/*Figure out size of buffer to be read: */
	fseek(restartfid, 0L, SEEK_END); 
	femmodel_size = ftell(restartfid);
	fseek(restartfid, 0L, SEEK_SET);

	/*Allocate buffer: */
	femmodel_buffer=xNew<char>(femmodel_size); 

	/*Read buffer from file: */
	fread_return=fread(femmodel_buffer,femmodel_size,sizeof(char),restartfid); if(fread_return!=1)_error_("error reading the buffer from marshalled file!");
	femmodel_buffer_ini=femmodel_buffer; //keep track of the initial position, so as to free later.

	/*Create new FemModel by demarshalling the buffer: */
	this->Marshall(&femmodel_buffer,NULL,MARSHALLING_BACKWARD);

	/*Reset position of buffer: */
	femmodel_buffer=femmodel_buffer_ini;

	/*Done, close file :*/
	pfclose(restartfid,restartfilename);

	/*Free ressources: */
	xDelete<char>(restartfilename);
	xDelete<char>(femmodel_buffer);
}/*}}}*/
void FemModel::SetCurrentConfiguration(int configuration_type,int analysis_type){/*{{{*/

	/*Use configuration_type to setup the analysis counter, the configurations of objects etc ... but use 
	 * analysis_type to drive the element numerics. This allows for use of 1 configuration_type for several 
	 * analyses. For example: do a SurfaceSlopeX, SurfaceSlopeY, BedSlopeX and BedSlopeY analysis using the 
	 * Slope configuration.*/
	int found=-1;
	for(int i=0;i<nummodels;i++){
	
		if (analysis_type_list[i]==configuration_type){
			found=i;
			break;
		}
	}
	if(found!=-1) analysis_counter=found;
	else _error_("Could not find alias for analysis_type " << EnumToStringx(configuration_type) << " in list of FemModel analyses");

	/*Now, plug analysis_counter and analysis_type inside the parameters: */
	this->parameters->SetParam(analysis_counter,AnalysisCounterEnum);
	this->parameters->SetParam(analysis_type,AnalysisTypeEnum);
	this->parameters->SetParam(configuration_type,ConfigurationTypeEnum);

	/*configure elements, loads and nodes, for this new analysis: */
	this->elements->SetCurrentConfiguration(elements,loads, nodes,vertices, materials,parameters);
	this->loads->SetCurrentConfiguration(elements, loads, nodes,vertices, materials,parameters);

	/*take care of toolkits options, that depend on this analysis type (present only after model processor)*/
	if(this->parameters->Exist(ToolkitsOptionsStringsEnum)){
		ToolkitsOptionsFromAnalysis(this->parameters,analysis_type);
		if(VerboseSolver()) _printf0_("      toolkits Options set for analysis type: " << EnumToStringx(analysis_type) << "\n");
	}

}
/*}}}*/
void FemModel::SetCurrentConfiguration(int configuration_type){/*{{{*/
	this->SetCurrentConfiguration(configuration_type,configuration_type);
}
/*}}}*/
int  FemModel::Size(){ /*{{{*/
	int   femmodel_size;

	this->Marshall(NULL,&femmodel_size,MARSHALLING_SIZE);

	return femmodel_size;
}
/*}}}*/
void FemModel::SolutionAnalysesList(int** panalyses,int* pnumanalyses,IoModel* iomodel,int solutiontype){/*{{{*/

	/*output: */
	int  numanalyses = 0;
	int* analyses=NULL;

	/*Intermediaries*/
	const int MAXANALYSES = 30;
	int   analyses_temp[MAXANALYSES];

	/*Analyses lists*/
	switch(solutiontype){

		case StressbalanceSolutionEnum:{
			bool isSIA,isFS;
			int  fe_FS;
			iomodel->FindConstant(&fe_FS,"md.flowequation.fe_FS");
			iomodel->FindConstant(&isSIA,"md.flowequation.isSIA");
			iomodel->FindConstant(&isFS,"md.flowequation.isFS");
			analyses_temp[numanalyses++]=StressbalanceAnalysisEnum;
			analyses_temp[numanalyses++]=StressbalanceVerticalAnalysisEnum;
			if(isSIA){
				analyses_temp[numanalyses++]=StressbalanceSIAAnalysisEnum;
			}
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			analyses_temp[numanalyses++]=ExtrudeFromBaseAnalysisEnum;
			analyses_temp[numanalyses++]=DepthAverageAnalysisEnum;
			if(fe_FS==LATaylorHoodEnum || fe_FS==LACrouzeixRaviartEnum){
				analyses_temp[numanalyses++]=UzawaPressureAnalysisEnum;
			}
			}
			break;

		case SteadystateSolutionEnum:{
			bool isSIA,isenthalpy;
			iomodel->FindConstant(&isSIA,"md.flowequation.isSIA");
			iomodel->FindConstant(&isenthalpy,"md.thermal.isenthalpy");
			analyses_temp[numanalyses++]=StressbalanceAnalysisEnum;
			analyses_temp[numanalyses++]=StressbalanceVerticalAnalysisEnum;
			if(isSIA){
				analyses_temp[numanalyses++]=StressbalanceSIAAnalysisEnum;
			}
			if(isenthalpy){
				analyses_temp[numanalyses++]=EnthalpyAnalysisEnum;
			}
			else{
				analyses_temp[numanalyses++]=ThermalAnalysisEnum;
				analyses_temp[numanalyses++]=MeltingAnalysisEnum;
			}
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			}
			break;

		case ThermalSolutionEnum:{
			bool isenthalpy;
			iomodel->FindConstant(&isenthalpy,"md.thermal.isenthalpy");
			if(isenthalpy){
				analyses_temp[numanalyses++]=EnthalpyAnalysisEnum;
			}
			else{
				analyses_temp[numanalyses++]=ThermalAnalysisEnum;
				analyses_temp[numanalyses++]=MeltingAnalysisEnum;
			}
			}
			break;

		case HydrologySolutionEnum:
			analyses_temp[numanalyses++]=HydrologyShreveAnalysisEnum;
			analyses_temp[numanalyses++]=HydrologyDCInefficientAnalysisEnum;
			analyses_temp[numanalyses++]=HydrologyDCEfficientAnalysisEnum;
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			analyses_temp[numanalyses++]=L2ProjectionEPLAnalysisEnum;
			break;

		case MasstransportSolutionEnum:
			analyses_temp[numanalyses++]=DepthAverageAnalysisEnum;
			analyses_temp[numanalyses++]=SmbAnalysisEnum;
			analyses_temp[numanalyses++]=MasstransportAnalysisEnum;
			analyses_temp[numanalyses++]=ExtrudeFromBaseAnalysisEnum;
			break;

		case BalancethicknessSolutionEnum:
			analyses_temp[numanalyses++]=BalancethicknessAnalysisEnum;
			break;

		case Balancethickness2SolutionEnum:
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			analyses_temp[numanalyses++]=SmoothAnalysisEnum;
			analyses_temp[numanalyses++]=Balancethickness2AnalysisEnum;
			break;

		case BalancethicknessSoftSolutionEnum:
			analyses_temp[numanalyses++]=BalancethicknessAnalysisEnum;
			break;

		case BalancevelocitySolutionEnum:
			analyses_temp[numanalyses++]=BalancevelocityAnalysisEnum;
			analyses_temp[numanalyses++]=SmoothAnalysisEnum;
			break;

		case SurfaceSlopeSolutionEnum:
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			break;

		case BedSlopeSolutionEnum:
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			break;

		case GiaSolutionEnum:
			analyses_temp[numanalyses++]=GiaIvinsAnalysisEnum;
			break;
		
		case EsaSolutionEnum:
			analyses_temp[numanalyses++]=EsaAnalysisEnum;
			break;
		
		case SealevelriseSolutionEnum:
			analyses_temp[numanalyses++]=SealevelriseAnalysisEnum;
			break;

		case SmbSolutionEnum:
			analyses_temp[numanalyses++]=SmbAnalysisEnum;
			break;

		case DamageEvolutionSolutionEnum:
			analyses_temp[numanalyses++]=DamageEvolutionAnalysisEnum;
			break;

		case TransientSolutionEnum:{
			bool isSIA,isFS,isthermal,isenthalpy,ismasstransport,isgroundingline,isstressbalance,ismovingfront,ishydrology,isdamage,issmb,isslr,isesa,isgia;
			iomodel->FindConstant(&isSIA,"md.flowequation.isSIA");
			iomodel->FindConstant(&isFS,"md.flowequation.isFS");
			iomodel->FindConstant(&isthermal,"md.transient.isthermal");
			iomodel->FindConstant(&isenthalpy,"md.thermal.isenthalpy");
			iomodel->FindConstant(&ismovingfront,"md.transient.ismovingfront");
			iomodel->FindConstant(&ismasstransport,"md.transient.ismasstransport");
			iomodel->FindConstant(&isstressbalance,"md.transient.isstressbalance");
			iomodel->FindConstant(&isgroundingline,"md.transient.isgroundingline");
			iomodel->FindConstant(&isdamage,"md.transient.isdamageevolution");
			iomodel->FindConstant(&ishydrology,"md.transient.ishydrology");
			iomodel->FindConstant(&issmb,"md.transient.issmb");
			iomodel->FindConstant(&isslr,"md.transient.isslr");
			iomodel->FindConstant(&isesa,"md.transient.isesa");
			iomodel->FindConstant(&isgia,"md.transient.isgia");
			if(isstressbalance){
				int  fe_FS;
				iomodel->FindConstant(&fe_FS,"md.flowequation.fe_FS");
				analyses_temp[numanalyses++]=StressbalanceAnalysisEnum;
				analyses_temp[numanalyses++]=StressbalanceVerticalAnalysisEnum;
				if(isSIA){
					analyses_temp[numanalyses++]=StressbalanceSIAAnalysisEnum;
				}
				analyses_temp[numanalyses++]=DepthAverageAnalysisEnum;
				if(fe_FS==LATaylorHoodEnum || fe_FS==LACrouzeixRaviartEnum){
					analyses_temp[numanalyses++]=UzawaPressureAnalysisEnum;
				}
			}
			if(isthermal && iomodel->domaintype==Domain3DEnum){
				if(isenthalpy){
					analyses_temp[numanalyses++]=EnthalpyAnalysisEnum;
				}
				else{
					analyses_temp[numanalyses++]=ThermalAnalysisEnum;
					analyses_temp[numanalyses++]=MeltingAnalysisEnum;
				}
			}
			if(ismasstransport || isgroundingline){
				analyses_temp[numanalyses++]=MasstransportAnalysisEnum;
			}
			if(issmb) analyses_temp[numanalyses++]=SmbAnalysisEnum;
			if(ismovingfront){
				analyses_temp[numanalyses++]=ExtrapolationAnalysisEnum;
				analyses_temp[numanalyses++]=LevelsetAnalysisEnum;
			}
			if(ishydrology){
				analyses_temp[numanalyses++]=HydrologyShreveAnalysisEnum;
				analyses_temp[numanalyses++]=HydrologySommersAnalysisEnum;
				analyses_temp[numanalyses++]=HydrologyDCInefficientAnalysisEnum;
				analyses_temp[numanalyses++]=HydrologyDCEfficientAnalysisEnum;
				analyses_temp[numanalyses++]=L2ProjectionEPLAnalysisEnum;
			}
			if(isdamage){
				analyses_temp[numanalyses++]=DamageEvolutionAnalysisEnum;
			}
			if(isslr){
				analyses_temp[numanalyses++]=SealevelriseAnalysisEnum;
			}
			if(isesa){
				analyses_temp[numanalyses++]=EsaAnalysisEnum;
			}
			if(isgia){
				analyses_temp[numanalyses++]=GiaIvinsAnalysisEnum;
			}

			if(iomodel->domaintype==Domain2DverticalEnum || iomodel->domaintype==Domain3DEnum){
				analyses_temp[numanalyses++]=ExtrudeFromBaseAnalysisEnum;
				analyses_temp[numanalyses++]=ExtrudeFromTopAnalysisEnum;
				analyses_temp[numanalyses++]=FreeSurfaceBaseAnalysisEnum;
				analyses_temp[numanalyses++]=FreeSurfaceTopAnalysisEnum;
			}
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			}
			break;

		default:
			_error_("solution type: " << EnumToStringx(solutiontype) << " not supported yet!");
			break;
	}

	/*Copy analyses from temp to output*/
	_assert_(numanalyses<MAXANALYSES);
	analyses=xNew<int>(numanalyses);
	for(int i=0;i<numanalyses;i++) analyses[i]=analyses_temp[i];

	/*Assign output pointers:*/
	if(pnumanalyses) *pnumanalyses=numanalyses;
	if(panalyses)    *panalyses=analyses;
	else              xDelete<int>(analyses);
}/*}}}*/
void FemModel::Solve(void){/*{{{*/

	/*profiling: */
	bool profiling = false;
	IssmDouble solution_time;
	IssmDouble solution_flops;
	IssmDouble solution_memory;

	/*solution: */
	int solution_type;
	void (*solutioncore)(FemModel*)=NULL; //core solution function pointer

	_printf0_("call computational core:\n");

	/*Retrieve solution_type from parameters: */
	parameters->FindParam(&solution_type,SolutionTypeEnum);

	/*Figure out which solution core we are going to run with the current solution type: */
	WrapperCorePointerFromSolutionEnum(&solutioncore,this->parameters,solution_type);

	/*run solution core: */
	profiler->Tag(STARTCORE);   
	solutioncore(this); 
	profiler->Tag(FINISHCORE);

	/*run AD core if needed: */
	profiler->Tag(STARTADCORE); 
	ad_core(this);      
	profiler->Tag(FINISHADCORE);

	/*some profiling results for the core: */
	parameters->FindParam(&profiling,DebugProfilingEnum);
	if(profiling){

		solution_time=profiler->DeltaTime(STARTCORE,FINISHCORE);
		solution_flops=profiler->DeltaFlops(STARTCORE,FINISHCORE);
		solution_memory=profiler->Memory(FINISHCORE);

		_printf0_("Core solution elapsed time    : " << solution_time   << " Seconds\n");
		_printf0_("Core solution number of flops : " << solution_flops  << " Flops\n");
		_printf0_("Core solution memory used     : " << solution_memory << " Bytes\n");

		/*Add to results: */
		results->AddObject(new GenericExternalResult<IssmDouble>(results->Size()+1, ProfilingSolutionTimeEnum,  solution_time));
		results->AddObject(new GenericExternalResult<IssmDouble>(results->Size()+1, ProfilingCurrentMemEnum,  solution_memory));
		results->AddObject(new GenericExternalResult<IssmDouble>(results->Size()+1, ProfilingCurrentFlopsEnum, solution_flops));

		#ifdef _HAVE_ADOLC_
		solution_time=profiler->DeltaTime(STARTADCORE,FINISHADCORE);
		solution_flops=profiler->DeltaFlops(STARTADCORE,FINISHADCORE);
		solution_memory=profiler->Memory(FINISHADCORE);

		_printf0_("AD Solution elapsed time    : " << solution_time   << " Seconds\n");
		_printf0_("AD Solution number of flops : " << solution_flops  << " Flops\n");
		_printf0_("AD Solution memory used     : " << solution_memory << " Bytes\n");
		#endif

	}
}
/*}}}*/

/*Modules:*/
void FemModel::BalancethicknessMisfitx(IssmDouble* presponse){/*{{{*/

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	IssmDouble  weight,vx,vy,H,dvx[2],dvy[2],dH[2];
	IssmDouble  temp,Jdet,dhdt,groundedice_melting,surface_mass_balance;
	IssmDouble* xyz_list = NULL;
	IssmDouble  dp[3];

	/*Compute Misfit: */
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));

		/*If on water, return 0: */
		if(!element->IsIceInElement()) continue;

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);
		Input* weights_input                   = element->GetInput(InversionCostFunctionsCoefficientsEnum);   _assert_(weights_input);
		Input* thickness_input                 = element->GetInput(ThicknessEnum); _assert_(thickness_input);
		Input* vx_input                        = element->GetInput(VxEnum);                                  _assert_(vx_input);
		Input* vy_input                        = element->GetInput(VyEnum);                                  _assert_(vy_input);
		Input* surface_mass_balance_input      = element->GetInput(SmbMassBalanceEnum);          _assert_(surface_mass_balance_input);
		Input* groundedice_melting_input       = element->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(groundedice_melting_input);
		Input* dhdt_input                      = element->GetInput(BalancethicknessThickeningRateEnum);      _assert_(dhdt_input);

		/* Start  looping on the number of gaussian points: */
		Gauss* gauss=element->NewGauss(2);
		for(int ig=gauss->begin();ig<gauss->end();ig++){

			gauss->GaussPoint(ig);

			/* Get Jacobian determinant: */
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			/*Get all parameters at gaussian point*/
			weights_input->GetInputValue(&weight,gauss,BalancethicknessMisfitEnum);
			thickness_input->GetInputValue(&H, gauss);
			thickness_input->GetInputDerivativeValue(&dH[0],xyz_list,gauss);
			surface_mass_balance_input->GetInputValue(&surface_mass_balance,gauss);
			groundedice_melting_input->GetInputValue(&groundedice_melting,gauss);
			dhdt_input->GetInputValue(&dhdt,gauss);
			vx_input->GetInputValue(&vx,gauss);
			vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
			vy_input->GetInputValue(&vy,gauss);
			vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);

			/*Balance thickness soft constraint J = 1/2 (div(Hv)-a)^2*/
			temp  = vx*dH[0]+vy*dH[1]+H*(dvx[0]+dvy[1]) - (surface_mass_balance-groundedice_melting-dhdt);
			J    +=weight*1/2*temp*temp*Jdet*gauss->weight;
		}

		/*clean up and Return: */
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*presponse=J;

}/*}}}*/
void FemModel::CalvingRateDevx(){/*{{{*/

	for(int i=0;i<elements->Size();i++){
		Element* element=dynamic_cast<Element*>(this->elements->GetObjectByOffset(i));
		element->CalvingRateDev();
	}
}
/*}}}*/
void FemModel::CalvingRateLevermannx(){/*{{{*/

	for(int i=0;i<elements->Size();i++){
		Element* element=dynamic_cast<Element*>(this->elements->GetObjectByOffset(i));
		element->CalvingRateLevermann();
	}
}
/*}}}*/
void FemModel::CostFunctionx(IssmDouble* pJ,IssmDouble** pJlist,int* pn){/*{{{*/

	/*Intermediary*/
	int      num_responses;
	int     *responses      = NULL;
	Results *cost_functions = NULL;

	/*Recover parameters*/
	parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	parameters->FindParam(&responses,NULL,InversionCostFunctionsEnum);

	/*Get the value of all cost functions*/
	this->RequestedOutputsx(&cost_functions,responses,num_responses);

	/*Get and add all contributions one by one*/
	IssmDouble  J=0;
	IssmDouble* Jlist = xNew<IssmDouble>(num_responses);
	for(int i=0;i<num_responses;i++){
		ExternalResult* result=(ExternalResult*)cost_functions->GetObjectByOffset(i);
		Jlist[i] = reCast<IssmDouble>(result->GetValue());
		J       += Jlist[i];
	}
	_assert_(cost_functions->Size()==num_responses);

	/*Assign output pointers: */
	delete cost_functions;
	xDelete<int>(responses);
	if(pJ)     *pJ     = J;
	if(pJlist) *pJlist = Jlist;
	else        xDelete<IssmDouble>(Jlist);
	if(pn)     *pn     = num_responses;
}
/*}}}*/
void FemModel::DeviatoricStressx(){/*{{{*/

	for(int i=0;i<elements->Size();i++){
		Element* element=dynamic_cast<Element*>(this->elements->GetObjectByOffset(i));
		element->ComputeDeviatoricStressTensor();
	}
}
/*}}}*/
void FemModel::Divergencex(IssmDouble* pdiv){/*{{{*/

	IssmDouble local_divergence=0;
	IssmDouble total_divergence;

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		local_divergence+=element->Divergence();
	}
	ISSM_MPI_Reduce(&local_divergence,&total_divergence,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_divergence,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pdiv=total_divergence;

}/*}}}*/
void FemModel::ElementOperationx(void (Element::*function)(void)){ /*{{{*/

	for(int i=0;i<elements->Size();i++){
		Element* element=dynamic_cast<Element*>(this->elements->GetObjectByOffset(i));
		(element->*function)();
	}

}
/*}}}*/
void FemModel::ElementResponsex(IssmDouble* presponse,int response_enum){/*{{{*/

	int found=0;
	int sumfound=0;
	int cpu_found=-1;
	int index;
	IssmDouble response;
	Element* element=NULL;

	/*retrieve element we are interested in: */
	this->parameters->FindParam(&index,IndexEnum);
	int my_rank=IssmComm::GetRank();

	/*now, go through our elements, and retrieve the one with this id: index: */
	for(int i=0;i<this->elements->Size();i++){
		element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		if (element->Id()==index){
			found=1;
			cpu_found=my_rank;
			break;
		}
	}

	/*Broadcast whether we found the element: */
	ISSM_MPI_Allreduce ( &found,&sumfound,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
	if(!sumfound)_error_("could not find material with id" << index << " to compute ElementResponse");

	/*Ok, we found the element, compute responseocity: */
	if(my_rank==cpu_found){
		element->ElementResponse(&response,response_enum);
	}

	/*Broadcast and plug into response: */
	ISSM_MPI_Allreduce ( &cpu_found,&cpu_found,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());
	ISSM_MPI_Bcast(&response,1,ISSM_MPI_DOUBLE,cpu_found,IssmComm::GetComm()); 

	/*Assign output pointers: */
	*presponse=response;

}/*}}}*/
void FemModel::FloatingAreax(IssmDouble* pV){/*{{{*/

	IssmDouble local_floating_area= 0;
	IssmDouble total_floating_area;

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		local_floating_area+=element->FloatingArea();
	}
	ISSM_MPI_Reduce(&local_floating_area,&total_floating_area,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_floating_area,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pV=total_floating_area;

}/*}}}*/
void FemModel::GetInputLocalMinMaxOnNodesx(IssmDouble** pmin,IssmDouble** pmax,IssmDouble* ug){/*{{{*/

	/*Get vector sizes for current configuration*/
	int configuration_type;
	this->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	int numnodes = this->nodes->NumberOfNodes(configuration_type);

	/*Initialize output vectors*/
	IssmDouble* uLmin_local = xNew<IssmDouble>(numnodes);
	IssmDouble* uLmax_local = xNew<IssmDouble>(numnodes);
	IssmDouble* uLmin = xNew<IssmDouble>(numnodes);
	IssmDouble* uLmax = xNew<IssmDouble>(numnodes);
	for(int i=0;i<numnodes;i++){
		uLmin_local[i] = +1.e+50;
		uLmax_local[i] = -1.e+50;
	}

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element->GetInputLocalMinMaxOnNodes(uLmin_local,uLmax_local,ug);
	}

	/*Synchronize all CPUs*/
	ISSM_MPI_Allreduce((void*)uLmin_local,(void*)uLmin,numnodes,ISSM_MPI_DOUBLE,ISSM_MPI_MIN,IssmComm::GetComm());
	ISSM_MPI_Allreduce((void*)uLmax_local,(void*)uLmax,numnodes,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,IssmComm::GetComm());
	xDelete<IssmDouble>(uLmin_local);
	xDelete<IssmDouble>(uLmax_local);

	/*Assign output pointers: */
	*pmin=uLmin;
	*pmax=uLmax;

}/*}}}*/
void FemModel::GroundedAreax(IssmDouble* pV){/*{{{*/

	IssmDouble local_grounded_area= 0;
	IssmDouble total_grounded_area;

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		local_grounded_area+=element->GroundedArea();
	}
	ISSM_MPI_Reduce(&local_grounded_area,&total_grounded_area,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_grounded_area,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pV=total_grounded_area;

}/*}}}*/
void FemModel::IceMassx(IssmDouble* pM){/*{{{*/

	IssmDouble local_ice_mass = 0;
	IssmDouble total_ice_mass;

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		local_ice_mass+=element->IceMass();
	}
	ISSM_MPI_Reduce(&local_ice_mass,&total_ice_mass,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_ice_mass,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pM=total_ice_mass;

}/*}}}*/
void FemModel::IceVolumeAboveFloatationx(IssmDouble* pV){/*{{{*/

	IssmDouble local_ice_volume_af = 0;
	IssmDouble total_ice_volume_af;

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		local_ice_volume_af+=element->IceVolumeAboveFloatation();
	}
	ISSM_MPI_Reduce(&local_ice_volume_af,&total_ice_volume_af,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_ice_volume_af,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pV=total_ice_volume_af;

}/*}}}*/
void FemModel::IceVolumex(IssmDouble* pV){/*{{{*/

	IssmDouble local_ice_volume = 0;
	IssmDouble total_ice_volume;

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		local_ice_volume+=element->IceVolume();
	}
	ISSM_MPI_Reduce(&local_ice_volume,&total_ice_volume,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_ice_volume,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pV=total_ice_volume;

}/*}}}*/
void FemModel::MassFluxx(IssmDouble* pmass_flux){/*{{{*/

	int          i,j;
	Element     *element       = NULL;
	int          element_id;
	bool         ispresent     = false;
	IssmDouble   mass_flux     = 0;
	IssmDouble   all_mass_flux = 0;
	int          counter;
	IssmDouble **array         = NULL;
	int          M;
	int         *mdims_array   = NULL;
	int         *ndims_array   = NULL;
	IssmDouble  *segments      = NULL;
	int          num_segments;

	/*First, figure out which segment to compute our mass flux on. Start with retrieving qmu_mass_flux_segments: */
	this->parameters->FindParam(&ispresent,MassFluxSegmentsPresentEnum);
	if(!ispresent)_error_("no mass flux segments available!");
	this->parameters->FindParam(&array,&M,&mdims_array,&ndims_array,MassFluxSegmentsEnum);

	/*Retrieve index of segments being used for MassFlux computation: */
	parameters->FindParam(&counter,IndexEnum);

	/*retrieve segments from array: */
	segments     = array[counter-1]; //matlab to "C" indexing
	num_segments = mdims_array[counter-1];

	/*Go through segments, and then elements, and figure out which elements belong to a segment. 
	 * When we find one, use the element to compute the mass flux on the segment: */
	for(i=0;i<num_segments;i++){
		element_id=reCast<int,IssmDouble>(*(segments+5*i+4));
		for(j=0;j<elements->Size();j++){
			element=(Element*)this->elements->GetObjectByOffset(j);
			if (element->Id()==element_id){
				/*We found the element which owns this segment, use it to compute the mass flux: */
				mass_flux+=element->MassFlux(segments+5*i+0);
				break;
			}
		}
	}

	ISSM_MPI_Allreduce ( (void*)&mass_flux,(void*)&all_mass_flux,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
	mass_flux=all_mass_flux;

	/*Free ressources:*/
	for(i=0;i<M;i++){
		IssmDouble* matrix=array[i];
		xDelete<IssmDouble>(matrix);
	}
	xDelete<int>(mdims_array);
	xDelete<int>(ndims_array);
	xDelete<IssmDouble*>(array);

	/*Assign output pointers: */
	*pmass_flux=mass_flux;

}/*}}}*/
void FemModel::MaxAbsVxx(IssmDouble* pmaxabsvx){/*{{{*/

	int i;
	IssmDouble maxabsvx;
	IssmDouble node_maxabsvx;
	IssmDouble element_maxabsvx;

	/*Go through elements, and request velocity: */
	maxabsvx=-INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_maxabsvx=element->inputs->MaxAbs(VxEnum);
		if(element_maxabsvx>maxabsvx) maxabsvx=element_maxabsvx;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxabsvx,&node_maxabsvx,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxabsvx,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	maxabsvx=node_maxabsvx;

	/*Assign output pointers:*/
	*pmaxabsvx=maxabsvx;

}/*}}}*/
void FemModel::MaxAbsVyx(IssmDouble* pmaxabsvy){/*{{{*/

	int i;
	IssmDouble maxabsvy;
	IssmDouble node_maxabsvy;
	IssmDouble element_maxabsvy;

	/*Go through elements, and request velocity: */
	maxabsvy=-INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_maxabsvy=element->inputs->MaxAbs(VyEnum);
		if(element_maxabsvy>maxabsvy) maxabsvy=element_maxabsvy;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxabsvy,&node_maxabsvy,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxabsvy,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	maxabsvy=node_maxabsvy;

	/*Assign output pointers:*/
	*pmaxabsvy=maxabsvy;

}/*}}}*/
void FemModel::MaxAbsVzx(IssmDouble* pmaxabsvz){/*{{{*/

	int i;
	IssmDouble maxabsvz;
	IssmDouble node_maxabsvz;
	IssmDouble element_maxabsvz;

	/*Go through elements, and request velocity: */
	maxabsvz=-INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_maxabsvz=element->inputs->MaxAbs(VzEnum);
		if(element_maxabsvz>maxabsvz) maxabsvz=element_maxabsvz;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxabsvz,&node_maxabsvz,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxabsvz,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	maxabsvz=node_maxabsvz;

	/*Assign output pointers:*/
	*pmaxabsvz=maxabsvz;

}/*}}}*/
void FemModel::MaxDivergencex(IssmDouble* pdiv){/*{{{*/

	IssmDouble local_divergence;
	IssmDouble node_max_divergence;
	IssmDouble max_divergence = -INFINITY;

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		local_divergence=element->Divergence();
		if(fabs(local_divergence)>max_divergence) max_divergence=fabs(local_divergence);
	}
	ISSM_MPI_Reduce(&max_divergence,&node_max_divergence,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_max_divergence,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	max_divergence=node_max_divergence;

	/*Assign output pointers: */
	*pdiv=max_divergence;

}/*}}}*/
void FemModel::MaxVelx(IssmDouble* pmaxvel){/*{{{*/

	int i;
	IssmDouble maxvel;
	IssmDouble node_maxvel;
	IssmDouble element_maxvel;

	/*Go through elements, and request velocity: */
	maxvel=-INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_maxvel = element->inputs->Max(VelEnum);
		if(element_maxvel>maxvel) maxvel=element_maxvel;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxvel,&node_maxvel,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxvel,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	maxvel=node_maxvel;

	/*Assign output pointers:*/
	*pmaxvel=maxvel;

}/*}}}*/
void FemModel::MaxVxx(IssmDouble* pmaxvx){/*{{{*/

	int i;
	IssmDouble maxvx;
	IssmDouble node_maxvx;
	IssmDouble element_maxvx;

	/*Go through elements, and request velocity: */
	maxvx=-INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_maxvx = element->inputs->Max(VxEnum);
		if(element_maxvx>maxvx) maxvx=element_maxvx;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxvx,&node_maxvx,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxvx,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	maxvx=node_maxvx;

	/*Assign output pointers:*/
	*pmaxvx=maxvx;

}/*}}}*/
void FemModel::MaxVyx(IssmDouble* pmaxvy){/*{{{*/

	int i;
	IssmDouble maxvy;
	IssmDouble node_maxvy;
	IssmDouble element_maxvy;

	/*Go through elements, and request velocity: */
	maxvy=-INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_maxvy = element->inputs->Max(VyEnum);
		if(element_maxvy>maxvy) maxvy=element_maxvy;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxvy,&node_maxvy,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxvy,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	maxvy=node_maxvy;

	/*Assign output pointers:*/
	*pmaxvy=maxvy;

}/*}}}*/
void FemModel::MaxVzx(IssmDouble* pmaxvz){/*{{{*/

	int i;
	IssmDouble maxvz;
	IssmDouble node_maxvz;
	IssmDouble element_maxvz;

	/*Go through elements, and request velocity: */
	maxvz=-INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_maxvz = element->inputs->Max(VzEnum);
		if(element_maxvz>maxvz) maxvz=element_maxvz;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxvz,&node_maxvz,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxvz,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	maxvz=node_maxvz;

	/*Assign output pointers:*/
	*pmaxvz=maxvz;

}/*}}}*/
void FemModel::MinVelx(IssmDouble* pminvel){/*{{{*/

	int i;
	IssmDouble minvel;
	IssmDouble node_minvel;
	IssmDouble element_minvel;

	/*Go through elements, and request velocity: */
	minvel=INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_minvel = element->inputs->Min(VelEnum);
		if(element_minvel<minvel) minvel=element_minvel;
	}

	/*Figure out minimum across the cluster: */
	ISSM_MPI_Reduce(&minvel,&node_minvel,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_minvel,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	minvel=node_minvel;

	/*Assign output pointers:*/
	*pminvel=minvel;

}/*}}}*/
void FemModel::MinVxx(IssmDouble* pminvx){/*{{{*/

	int i;
	IssmDouble minvx;
	IssmDouble node_minvx;
	IssmDouble element_minvx;

	/*Go through elements, and request velocity: */
	minvx=INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_minvx = element->inputs->Min(VxEnum);
		if(element_minvx<minvx) minvx=element_minvx;
	}

	/*Figure out minimum across the cluster: */
	ISSM_MPI_Reduce(&minvx,&node_minvx,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_minvx,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	minvx=node_minvx;

	/*Assign output pointers:*/
	*pminvx=minvx;

}/*}}}*/
void FemModel::MinVyx(IssmDouble* pminvy){/*{{{*/

	int i;
	IssmDouble minvy;
	IssmDouble node_minvy;
	IssmDouble element_minvy;

	/*Go through elements, and request velocity: */
	minvy=INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_minvy = element->inputs->Min(VyEnum);
		if(element_minvy<minvy) minvy=element_minvy;
	}

	/*Figure out minimum across the cluster: */
	ISSM_MPI_Reduce(&minvy,&node_minvy,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_minvy,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	minvy=node_minvy;

	/*Assign output pointers:*/
	*pminvy=minvy;

}/*}}}*/
void FemModel::MinVzx(IssmDouble* pminvz){/*{{{*/

	int i;
	IssmDouble minvz;
	IssmDouble node_minvz;
	IssmDouble element_minvz;

	/*Go through elements, and request velocity: */
	minvz=INFINITY;
	for(i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element_minvz = element->inputs->Min(VzEnum);
		if(element_minvz<minvz) minvz=element_minvz;
	}

	/*Figure out minimum across the cluster: */
	ISSM_MPI_Reduce(&minvz,&node_minvz,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_minvz,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());   
	minvz=node_minvz;

	/*Assign output pointers:*/
	*pminvz=minvz;

}/*}}}*/
void FemModel::OutputControlsx(Results **presults){/*{{{*/

	/*parameters: */
	int         num_controls,step;
	IssmDouble  time;
	int        *control_type = NULL;

	/*recover results*/
	Results* results = *presults;
	if(!results) results = new Results();

	/*Get list of Controls*/
	this->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	this->parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);
	this->parameters->FindParam(&step,StepEnum);
	this->parameters->FindParam(&time,TimeEnum);

	for(int i=0;i<num_controls;i++){

		int control_enum = control_type[i];
		int gradient_enum;

		switch(i){
			case 0: gradient_enum = Gradient1Enum; break;
			case 1: gradient_enum = Gradient2Enum; break;
			case 2: gradient_enum = Gradient3Enum; break;
			default: _error_("more than 3 controls not implemented yet");
		}

		/*Allocate vector*/
		Vector<IssmPDouble> *vector_control  = new Vector<IssmPDouble>(this->vertices->NumberOfVertices());
		Vector<IssmPDouble> *vector_gradient = new Vector<IssmPDouble>(this->vertices->NumberOfVertices());

		/*Fill in vector*/
		for(int j=0;j<elements->Size();j++){
			Element* element=(Element*)elements->GetObjectByOffset(j);
			element->ControlToVectors(vector_control,vector_gradient,control_enum);
		}
		vector_control->Assemble();
		vector_gradient->Assemble();

		results->AddResult(new GenericExternalResult<Vector<IssmPDouble>*>(results->Size()+1,control_enum,vector_control ,step,time));
		results->AddResult(new GenericExternalResult<Vector<IssmPDouble>*>(results->Size()+1,gradient_enum,vector_gradient,step,time));
	}

	/*Clean up and return*/
	xDelete<int>(control_type);
}
/*}}}*/
void FemModel::RequestedDependentsx(void){/*{{{*/

	bool        isautodiff      = false;
	IssmDouble  output_value;

	int         num_dependents;
	IssmPDouble *dependents;
	DataSet*    dependent_objects=NULL;
	int my_rank=IssmComm::GetRank();

	/*AD mode on?: */
	parameters->FindParam(&isautodiff,AutodiffIsautodiffEnum);

	if(isautodiff){
		#ifdef _HAVE_ADOLC_
		parameters->FindParam(&num_dependents,AutodiffNumDependentsEnum);
		parameters->FindParam(&dependent_objects,AutodiffDependentObjectsEnum);
		if(num_dependents){
			dependents=xNew<IssmPDouble>(num_dependents);

			/*Go through our dependent variables, and compute the response:*/
			for(int i=0;i<dependent_objects->Size();i++){
				DependentObject* dep=(DependentObject*)dependent_objects->GetObjectByOffset(i);
				dep->Responsex(&output_value,this);
				if (my_rank==0) {
					output_value>>=dependents[i];
				}
			}
		}
		delete dependent_objects;
		if(num_dependents)xDelete<IssmPDouble>(dependents);
		#else
		_error_("Should not be requesting dependents when an AD library is not available!");
		#endif
	}
}
/*}}}*/
void FemModel::RequestedOutputsx(Results **presults,char** requested_outputs, int numoutputs, bool save_results){/*{{{*/

	/*Intermediaries*/
	bool        isvec,results_on_nodes;
	int         step,output_enum;
	IssmDouble  time;
	IssmDouble  double_result;
	const char *output_string = NULL;

	/*recover results*/
	Results* results = *presults;
	if(!results) results = new Results();

	/*Get time and step*/
	parameters->FindParam(&step,StepEnum);
	parameters->FindParam(&time,TimeEnum);
	parameters->FindParam(&results_on_nodes,SettingsResultsOnNodesEnum);

	/*Go through all requested output*/
	for(int i=0;i<numoutputs;i++){

		output_string = requested_outputs[i];
		output_enum   = StringToEnumx(output_string,false);
		isvec         = false;

		/*If string is not an enum, it is defined in output definitions*/
		if(output_enum<0){
			double_result = OutputDefinitionsResponsex(this,output_string);
			if(save_results){
				results->AddResult(new GenericExternalResult<IssmPDouble>(results->Size()+1,output_string,reCast<IssmPDouble>(double_result),step,time));
				continue;
			}
		}
		else{
			/*last chance for the output definition, if the enum is one of Outputdefinition[1-10]Enum:*/
			if(output_enum>=Outputdefinition1Enum && output_enum <=Outputdefinition10Enum){
				double_result = OutputDefinitionsResponsex(this,output_enum);
				if(save_results){
					results->AddResult(new GenericExternalResult<IssmPDouble>(results->Size()+1,output_string,reCast<IssmPDouble>(double_result),step,time));
					continue;
				}
			}
			else{
				switch(output_enum){

					/*Scalar output*/
					case DivergenceEnum:               this->Divergencex(&double_result);               break;
					case MaxDivergenceEnum:            this->MaxDivergencex(&double_result);            break;
					case IceMassEnum:                  this->IceMassx(&double_result);                  break;
					case IceVolumeEnum:                this->IceVolumex(&double_result);                break;
					case IceVolumeAboveFloatationEnum: this->IceVolumeAboveFloatationx(&double_result); break;
					case GroundedAreaEnum:             this->GroundedAreax(&double_result);             break;
					case FloatingAreaEnum:             this->FloatingAreax(&double_result);             break;
					case MinVelEnum:                   this->MinVelx(&double_result);                   break;
					case MaxVelEnum:                   this->MaxVelx(&double_result);                   break;
					case MinVxEnum:                    this->MinVxx(&double_result);                    break;
					case MaxVxEnum:                    this->MaxVxx(&double_result);                    break;
					case MaxAbsVxEnum:                 this->MaxAbsVxx(&double_result);                 break;
					case MinVyEnum:                    this->MinVyx(&double_result);                    break;
					case MaxVyEnum:                    this->MaxVyx(&double_result);                    break;
					case MaxAbsVyEnum:                 this->MaxAbsVyx(&double_result);                 break;
					case MinVzEnum:                    this->MinVzx(&double_result);                    break;
					case MaxVzEnum:                    this->MaxVzx(&double_result);                    break;
					case MaxAbsVzEnum:                 this->MaxAbsVzx(&double_result);                 break;
					case MassFluxEnum:                 this->MassFluxx(&double_result);                 break;
					case TotalFloatingBmbEnum:         this->TotalFloatingBmbx(&double_result);         break;
					case TotalGroundedBmbEnum:         this->TotalGroundedBmbx(&double_result);         break;
					case TotalSmbEnum:                 this->TotalSmbx(&double_result);                 break;

			   /*Scalar control output*/
				case SurfaceAbsVelMisfitEnum:       SurfaceAbsVelMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);        break;
				case SurfaceRelVelMisfitEnum:       SurfaceRelVelMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);        break;
				case SurfaceLogVelMisfitEnum:       SurfaceLogVelMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);        break;
				case SurfaceLogVxVyMisfitEnum:      SurfaceLogVxVyMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);       break;
				case SurfaceAverageVelMisfitEnum:   SurfaceAverageVelMisfitx(&double_result,this);                                                  break;
				case ThicknessAbsMisfitEnum:        ThicknessAbsMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);         break;
				case ThicknessAbsGradientEnum:      this->ThicknessAbsGradientx(&double_result);                                                    break;
				case ThicknessAlongGradientEnum:    ThicknessAlongGradientx(&double_result,elements,nodes,vertices,loads,materials,parameters);     break;
				case ThicknessAcrossGradientEnum:   ThicknessAcrossGradientx(&double_result,elements,nodes,vertices,loads,materials,parameters);    break;
				case ThicknessPositiveEnum:         this->ThicknessPositivex(&double_result);                                                       break;
				case RheologyBbarAbsGradientEnum:   RheologyBbarAbsGradientx(&double_result,elements,nodes,vertices,loads,materials,parameters);    break;
				case RheologyBAbsGradientEnum:      RheologyBAbsGradientx(&double_result,elements,nodes,vertices,loads,materials,parameters);       break;
				case DragCoefficientAbsGradientEnum:DragCoefficientAbsGradientx(&double_result,elements,nodes,vertices,loads,materials,parameters); break;
				case BalancethicknessMisfitEnum:    BalancethicknessMisfitx(&double_result);                                                        break;
				case SurfaceAbsMisfitEnum:          SurfaceAbsMisfitx(&double_result); break;

				   /*Vector */
					default:

						/*Vector layout*/
						int interpolation,nodesperelement,size,nlines,ncols,array_size;
						int rank_interpolation=-1,rank_nodesperelement=-1,rank_arraysize=-1,max_rank_arraysize=0;
						bool isarray=false;

						/*Get interpolation (and compute input if necessary)*/
						for(int j=0;j<elements->Size();j++){
							Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(j));
							element->ResultInterpolation(&rank_interpolation,&rank_nodesperelement,&rank_arraysize,output_enum);
							if (rank_arraysize>max_rank_arraysize)max_rank_arraysize=rank_arraysize;
						}
						rank_arraysize=max_rank_arraysize;

						/*Broadcast for cpus that do not have any elements*/
						ISSM_MPI_Reduce(&rank_interpolation,&interpolation,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
						ISSM_MPI_Reduce(&rank_nodesperelement,&nodesperelement,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
						ISSM_MPI_Reduce(&rank_arraysize,&array_size,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
						ISSM_MPI_Bcast(&interpolation,1,ISSM_MPI_INT,0,IssmComm::GetComm());
						ISSM_MPI_Bcast(&nodesperelement,1,ISSM_MPI_INT,0,IssmComm::GetComm());
						ISSM_MPI_Bcast(&array_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());
						
						if(results_on_nodes){

							/*Allocate matrices*/
							int         nbe       = this->elements->NumberOfElements();
							IssmDouble* values    = xNewZeroInit<IssmDouble>(nbe*nodesperelement);
							IssmDouble* allvalues = xNew<IssmDouble>(nbe*nodesperelement);

							/*Fill-in matrix*/
							for(int j=0;j<elements->Size();j++){
								Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(j));
								element->ResultToPatch(values,nodesperelement,output_enum);
							}

							/*Gather from all cpus*/
							ISSM_MPI_Allreduce((void*)values,(void*)allvalues,nbe*nodesperelement,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
							xDelete<IssmDouble>(values);

							if(save_results)results->AddResult(new GenericExternalResult<IssmDouble*>(results->Size()+1,output_enum,allvalues,nbe,nodesperelement,step,time));
							xDelete<IssmDouble>(allvalues);

						}
						else{

							/*Allocate vector depending on interpolation*/
							switch(interpolation){
								case P0Enum: isarray = false; size = this->elements->NumberOfElements(); break;
								case P1Enum: isarray = false; size = this->vertices->NumberOfVertices(); break;
								case P0ArrayEnum: isarray = true; nlines = this->elements->NumberOfElements(); ncols= array_size; break;
								default:     _error_("Interpolation "<<EnumToStringx(interpolation)<<" not supported yet");

							}
							if (!isarray){
								Vector<IssmDouble> *vector_result = new Vector<IssmDouble>(size);

								/*Fill in vector*/
								for(int j=0;j<elements->Size();j++){
									Element* element=(Element*)elements->GetObjectByOffset(j);
									element->ResultToVector(vector_result,output_enum);
								}
								vector_result->Assemble();

								if(save_results)results->AddResult(new GenericExternalResult<Vector<IssmDouble>*>(results->Size()+1,output_enum,vector_result,step,time));
							}
							else{
								IssmDouble* values    = xNewZeroInit<IssmDouble>(nlines*ncols);
								IssmDouble* allvalues = xNew<IssmDouble>(nlines*ncols);
								
								/*Fill-in matrix*/
								for(int j=0;j<elements->Size();j++){
									Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(j));
									element->ResultToMatrix(values,ncols, output_enum);
								}
								/*Gather from all cpus*/
								ISSM_MPI_Allreduce((void*)values,(void*)allvalues,ncols*nlines,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
								xDelete<IssmDouble>(values);
								
								if(save_results)results->AddResult(new GenericExternalResult<IssmDouble*>(results->Size()+1,output_enum,allvalues,nlines,ncols,step,time));
								xDelete<IssmDouble>(allvalues);
							}
						}
						isvec = true;
						break;
				}
			}

		}

		/*Add result to Results*/
		if(!isvec && save_results){
			results->AddResult(new GenericExternalResult<IssmPDouble>(results->Size()+1,output_string,reCast<IssmPDouble>(double_result),step,time));
		}
	}

	/*Assign pointer and clean up*/
	*presults = results;
}
/*}}}*/
void FemModel::RequestedOutputsx(Results **presults,int* requested_outputs, int numoutputs,bool save_results){/*{{{*/

	/*Convert list of enums to list of string*/
	char** enumlist = xNew<char*>(numoutputs);
	for(int i=0;i<numoutputs;i++){
		EnumToStringx(&enumlist[i],requested_outputs[i]);
	}

	/*Call main module*/
	this->RequestedOutputsx(presults,enumlist,numoutputs,save_results);

	/*clean up and return*/
	for(int i=0;i<numoutputs;i++) xDelete<char>(enumlist[i]);
	xDelete<char*>(enumlist);
	return;
}
/*}}}*/
void FemModel::ResetLevelset(void){/*{{{*/

	/*recover my_rank:*/
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*Get domain type (2d or 3d)*/
	int domaintype;
	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	
	/*1: go throug all elements of this partition and figure out how many
	 * segments we have (corresopnding to levelset = 0)*/
	DataSet* segments=new DataSet();
	for(int i=0;i<elements->Size();i++){
		Element* element=dynamic_cast<Element*>(this->elements->GetObjectByOffset(i));
		if(!element->IsOnBase()) continue;
		Element* basalelement = element->SpawnBasalElement();
		basalelement->WriteLevelsetSegment(segments);
		if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	}

	/*2: now get the segments from all partitions*/
	int  segcount=segments->Size();
	int* allsegcount=xNew<int>(num_procs);
	ISSM_MPI_Gather(&segcount,1,ISSM_MPI_INT,allsegcount,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(allsegcount,num_procs,ISSM_MPI_INT,0,IssmComm::GetComm());

	/* Every cpu should start its own dof count at the end of the dofcount from cpu-1*/
	int numseg_offset=0;
	int numseg=0;
	for(int i=0;i<my_rank;  i++) numseg_offset+=allsegcount[i];
	for(int i=0;i<num_procs;i++) numseg+=allsegcount[i];
	IssmDouble* segmentlist    = xNewZeroInit<IssmDouble>(4*numseg);
	IssmDouble* allsegmentlist = xNewZeroInit<IssmDouble>(4*numseg);
	for(int i=0;i<segments->Size();i++){
		Contour<IssmDouble>* segment=(Contour<IssmDouble>*)segments->GetObjectByOffset(i);
		_assert_(segment->nods == 2);
		segmentlist[(numseg_offset+i)*4 + 0] = segment->x[0];
		segmentlist[(numseg_offset+i)*4 + 1] = segment->y[0];
		segmentlist[(numseg_offset+i)*4 + 2] = segment->x[1];
		segmentlist[(numseg_offset+i)*4 + 3] = segment->y[1];
	}
	ISSM_MPI_Allreduce((void*)segmentlist,(void*)allsegmentlist,4*numseg,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
	delete segments;
	xDelete<IssmDouble>(segmentlist);
	xDelete<int>(allsegcount);

	/*3: update levelset for all elements*/
	for(int i=0;i<elements->Size();i++){
		Element* element=dynamic_cast<Element*>(this->elements->GetObjectByOffset(i));
		if(!element->IsOnBase()) continue;
		element->ResetLevelsetFromSegmentlist(allsegmentlist,numseg);
	}


	/*Extrude if necessary*/
	int elementtype;
	this->parameters->FindParam(&elementtype,MeshElementtypeEnum);
	if(elementtype==PentaEnum){
		InputExtrudex(this,MaskIceLevelsetEnum,-1);
	}
	else if(elementtype==TriaEnum){
		/*no need to extrude*/
	}
	else{
		_error_("not implemented yet");
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(allsegmentlist);

}
/*}}}*/
void FemModel::Responsex(IssmDouble* responses,const char* response_descriptor){/*{{{*/

	int response_descriptor_enum=StringToEnumx(response_descriptor);
	this->Responsex(responses, response_descriptor_enum);

}
/*}}}*/
void FemModel::Responsex(IssmDouble* responses,int response_descriptor_enum){/*{{{*/

	switch (response_descriptor_enum){

		case DivergenceEnum:               this->Divergencex(responses); break;
		case MaxDivergenceEnum:            this->MaxDivergencex(responses); break;
		case IceMassEnum:                  this->IceMassx(responses); break;
		case IceVolumeEnum:                this->IceVolumex(responses); break;
		case IceVolumeAboveFloatationEnum: this->IceVolumeAboveFloatationx(responses); break;
		case GroundedAreaEnum:             this->GroundedAreax(responses); break;
		case FloatingAreaEnum:             this->FloatingAreax(responses); break;
		case MinVelEnum:                   this->MinVelx(responses); break;
		case MaxVelEnum:                   this->MaxVelx(responses); break;
		case MinVxEnum:                    this->MinVxx(responses); break;
		case MaxVxEnum:                    this->MaxVxx(responses); break;
		case MaxAbsVxEnum:                 this->MaxAbsVxx(responses); break;
		case MinVyEnum:                    this->MinVyx(responses); break;
		case MaxVyEnum:                    this->MaxVyx(responses); break;
		case MaxAbsVyEnum:                 this->MaxAbsVyx(responses); break;
		case MinVzEnum:                    this->MinVzx(responses); break;
		case MaxVzEnum:                    this->MaxVzx(responses); break;
		case MaxAbsVzEnum:                 this->MaxAbsVzx(responses); break;
		case MassFluxEnum:                 this->MassFluxx(responses); break;
		case SurfaceAbsVelMisfitEnum:      SurfaceAbsVelMisfitx(responses, elements,nodes, vertices, loads, materials,parameters); break;
		case SurfaceRelVelMisfitEnum:      SurfaceRelVelMisfitx(responses, elements,nodes, vertices, loads, materials,parameters); break;
		case SurfaceLogVelMisfitEnum:      SurfaceLogVelMisfitx(responses, elements,nodes, vertices, loads, materials,parameters); break;
		case SurfaceLogVxVyMisfitEnum:     SurfaceLogVxVyMisfitx(responses, elements,nodes, vertices, loads, materials,parameters); break;
		case SurfaceAverageVelMisfitEnum:  SurfaceAverageVelMisfitx(responses,this); break;
		case ThicknessAbsMisfitEnum:       ThicknessAbsMisfitx(responses, elements,nodes, vertices, loads, materials, parameters); break;
		case ThicknessAbsGradientEnum:     this->ThicknessAbsGradientx(responses); break;
		case ThicknessAlongGradientEnum:   ThicknessAlongGradientx(responses, elements,nodes, vertices, loads, materials, parameters); break;
		case ThicknessAcrossGradientEnum:  ThicknessAcrossGradientx(responses, elements,nodes, vertices, loads, materials, parameters); break;
		case RheologyBbarAbsGradientEnum:  RheologyBbarAbsGradientx(responses, elements,nodes, vertices, loads, materials, parameters); break;
		case DragCoefficientAbsGradientEnum:DragCoefficientAbsGradientx(responses, elements,nodes, vertices, loads, materials, parameters); break;
		case BalancethicknessMisfitEnum:   BalancethicknessMisfitx(responses); break;
		case TotalFloatingBmbEnum:			  this->TotalFloatingBmbx(responses); break;
		case TotalGroundedBmbEnum:			  this->TotalGroundedBmbx(responses); break;
		case TotalSmbEnum:					  this->TotalSmbx(responses); break;
		case MaterialsRheologyBbarEnum:    this->ElementResponsex(responses,MaterialsRheologyBbarEnum); break;
		case VelEnum:                      this->ElementResponsex(responses,VelEnum); break;
		case FrictionCoefficientEnum:      NodalValuex(responses, FrictionCoefficientEnum,elements,nodes, vertices, loads, materials, parameters); break;
		default: _error_("response descriptor \"" << EnumToStringx(response_descriptor_enum) << "\" not supported yet!"); break; 
	}

}
/*}}}*/
void FemModel::StrainRateparallelx(){/*{{{*/

	for(int i=0;i<elements->Size();i++){
		Element* element=dynamic_cast<Element*>(this->elements->GetObjectByOffset(i));
		element->StrainRateparallel();
	}
}
/*}}}*/
void FemModel::StrainRateperpendicularx(){/*{{{*/

	for(int i=0;i<elements->Size();i++){
		Element* element=dynamic_cast<Element*>(this->elements->GetObjectByOffset(i));
		element->StrainRateperpendicular();
	}
}
/*}}}*/
void FemModel::StressIntensityFactorx(){/*{{{*/

	/*Update input for basal element only*/
	for(int i=0;i<elements->Size();i++){
		Element* element=dynamic_cast<Element*>(this->elements->GetObjectByOffset(i));
		element->StressIntensityFactor();
	}
}
	/*}}}*/
void FemModel::SurfaceAbsMisfitx(IssmDouble* presponse){/*{{{*/

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	IssmDouble  surface,surfaceobs,weight;
	IssmDouble  Jdet;
	IssmDouble* xyz_list = NULL;

	/*Compute Misfit: */
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));

		/*If on water, return 0: */
		if(!element->IsIceInElement()) continue;

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);

		 /*Retrieve all inputs we will be needing: */
		 Input* weights_input   =element->GetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
		 Input* surface_input   =element->GetInput(SurfaceEnum);                            _assert_(surface_input);
		 Input* surfaceobs_input=element->GetInput(InversionSurfaceObsEnum);                _assert_(surfaceobs_input);

		 /* Start  looping on the number of gaussian points: */
		 Gauss* gauss=element->NewGauss(2);
		 for(int ig=gauss->begin();ig<gauss->end();ig++){

			 gauss->GaussPoint(ig);

			 /* Get Jacobian determinant: */
			 element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			 /*Get all parameters at gaussian point*/
			 weights_input->GetInputValue(&weight,gauss,SurfaceAbsMisfitEnum);
			 surface_input->GetInputValue(&surface,gauss);
			 surfaceobs_input->GetInputValue(&surfaceobs,gauss);

			 /*Compute SurfaceAbsMisfitEnum*/
			 J+=0.5*(surface-surfaceobs)*(surface-surfaceobs)*weight*Jdet*gauss->weight;
		 }
		 delete gauss;
		 xDelete<IssmDouble>(xyz_list);
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*presponse=J;

}/*}}}*/
void FemModel::ThicknessAbsGradientx( IssmDouble* pJ){/*{{{*/

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	IssmDouble  thickness,thicknessobs,weight;
	IssmDouble  Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble  dp[3];

	/*Compute Misfit: */
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));

		/*If on water, return 0: */
		if(!element->IsIceInElement()) continue;

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);

		/*Retrieve all inputs we will be needing: */
		Input* weights_input   =element->GetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
		Input* thickness_input =element->GetInput(ThicknessEnum);                          _assert_(thickness_input);

		/* Start  looping on the number of gaussian points: */
		Gauss* gauss=element->NewGauss(2);
		for(int ig=gauss->begin();ig<gauss->end();ig++){

			gauss->GaussPoint(ig);

			/* Get Jacobian determinant: */
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			/*Get all parameters at gaussian point*/
			weights_input->GetInputValue(&weight,gauss,ThicknessAbsGradientEnum);
			thickness_input->GetInputDerivativeValue(&dp[0],xyz_list,gauss);

			/*Tikhonov regularization: J = 1/2 ((dp/dx)^2 + (dp/dy)^2) */ 
			J+=weight*1/2*(dp[0]*dp[0]+dp[1]*dp[1])*Jdet*gauss->weight;
		}

		/*clean up and Return: */
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}
/*}}}*/
void FemModel::ThicknessPositivex(IssmDouble* pJ){/*{{{*/

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	IssmDouble  thickness,weight;
	IssmDouble  Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble  H;

	/*Compute Misfit: */
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));

		/*If on water, return 0: */
		if(!element->IsIceInElement()) continue;

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);

		/*Retrieve all inputs we will be needing: */
		Input* weights_input   =element->GetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
		Input* thickness_input =element->GetInput(ThicknessEnum);                          _assert_(thickness_input);

		/* Start  looping on the number of gaussian points: */
		Gauss* gauss=element->NewGauss(2);
		for(int ig=gauss->begin();ig<gauss->end();ig++){

			gauss->GaussPoint(ig);

			/* Get Jacobian determinant: */
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			/*Get all parameters at gaussian point*/
			weights_input->GetInputValue(&weight,gauss,ThicknessPositiveEnum);
			thickness_input->GetInputValue(&H,gauss);

			/*int min(H,0)^2 */
			if(H<=0){
				J+=weight*H*H*Jdet*gauss->weight;
			}
		}

		/*clean up and Return: */
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}
/*}}}*/
void FemModel::TimeAdaptx(IssmDouble* pdt){/*{{{*/

	int      i;

	/*output: */
	IssmDouble   dt;

	/*intermediary: */
	Element *element     = NULL;
	IssmDouble   min_dt      = 0;
	IssmDouble   node_min_dt = 0;

	/*Go through elements, and figure out the minimum of the time steps for each element (using CFL criterion): */
	element=(Element*)elements->GetObjectByOffset(0); min_dt=element->TimeAdapt();

	for (i=1;i<elements->Size();i++){
		element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		dt=element->TimeAdapt();
		if(dt<min_dt)min_dt=dt;
	}

	/*Figure out minimum across the cluster: */
	ISSM_MPI_Reduce (&min_dt,&node_min_dt,1,ISSM_MPI_DOUBLE,ISSM_MPI_MIN,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_min_dt,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	min_dt=node_min_dt;

	/*Assign output pointers:*/
	*pdt=min_dt;
}
/*}}}*/
void FemModel::TotalFloatingBmbx(IssmDouble* pFbmb){/*{{{*/

	IssmDouble local_fbmb = 0;
	IssmDouble total_fbmb;

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		local_fbmb+=element->TotalFloatingBmb();
	}
	ISSM_MPI_Reduce(&local_fbmb,&total_fbmb,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_fbmb,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pFbmb=total_fbmb;

}/*}}}*/
void FemModel::TotalGroundedBmbx(IssmDouble* pGbmb){/*{{{*/

	IssmDouble local_gbmb = 0;
	IssmDouble total_gbmb;

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		local_gbmb+=element->TotalGroundedBmb();
	}
	ISSM_MPI_Reduce(&local_gbmb,&total_gbmb,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_gbmb,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pGbmb=total_gbmb;

}/*}}}*/
void FemModel::TotalSmbx(IssmDouble* pSmb){/*{{{*/

	IssmDouble local_smb = 0;
	IssmDouble total_smb;

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		local_smb+=element->TotalSmb();
	}
	ISSM_MPI_Reduce(&local_smb,&total_smb,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_smb,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pSmb=total_smb;

}/*}}}*/
void FemModel::UpdateConstraintsExtrudeFromBasex(void){ /*{{{*/

	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->UpdateConstraintsExtrudeFromBase();
	}

}
/*}}}*/
void FemModel::UpdateConstraintsExtrudeFromTopx(void){ /*{{{*/

	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->UpdateConstraintsExtrudeFromTop();
	}

}
/*}}}*/
void FemModel::UpdateConstraintsx(void){ /*{{{*/

	IssmDouble time;
	int        analysis_type;

	/*retrieve parameters: */
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	parameters->FindParam(&time,TimeEnum);

	/*start module: */
	if(VerboseModule()) _printf0_("   Updating constraints and active domain of analysis " << EnumToStringx(analysis_type)  << " for time: " << time << "\n");

	Analysis* analysis= EnumToAnalysis(analysis_type);
	analysis->UpdateConstraints(this);
	delete analysis;
	
	/*Second, constraints might be time dependent: */
	SpcNodesx(nodes,constraints,parameters,analysis_type); 

	/*Now, update degrees of freedoms: */
	NodesDofx(nodes,parameters,analysis_type);

}
/*}}}*/
int  FemModel::UpdateVertexPositionsx(void){ /*{{{*/

	IssmDouble         *surface = NULL;
	IssmDouble         *bed     = NULL;
			
	if(VerboseSolution()) _printf0_("   updating vertices positions\n");

	/*get vertex vectors for bed and thickness: */
	GetVectorFromInputsx(&surface  ,this, SurfaceEnum,VertexPIdEnum);
	GetVectorFromInputsx(&bed      ,this, BaseEnum,   VertexPIdEnum);

	/*Allocate vector*/
	Vector<IssmDouble> *vx=new Vector<IssmDouble>(vertices->NumberOfVertices());
	Vector<IssmDouble> *vy=new Vector<IssmDouble>(vertices->NumberOfVertices());
	Vector<IssmDouble> *vz=new Vector<IssmDouble>(vertices->NumberOfVertices());

	/*Update verices new geometry: */
	for(int i=0;i<vertices->Size();i++){
		Vertex* vertex=(Vertex*)vertices->GetObjectByOffset(i);
		vertex->UpdatePosition(vx,vy,vz,parameters,surface,bed);
	}

	/*Assemble mesh velocity*/
	vx->Assemble();
	vy->Assemble();
	vz->Assemble();

	/*Update element inputs*/
	InputUpdateFromVectorx(this,vx,VxMeshEnum,VertexPIdEnum);
	InputUpdateFromVectorx(this,vy,VyMeshEnum,VertexPIdEnum);
	InputUpdateFromVectorx(this,vz,VzMeshEnum,VertexPIdEnum);

	/*Free ressources:*/
	delete vx;
	delete vy;
	delete vz;
	xDelete<IssmDouble>(bed);
	xDelete<IssmDouble>(surface);
	return 1;
}
/*}}}*/
#ifdef  _HAVE_DAKOTA_
void FemModel::DakotaResponsesx(double* d_responses,char** responses_descriptors,int numresponsedescriptors,int d_numresponses){/*{{{*/

	int        i,j;
	int        my_rank;

	/*intermediary: */
	char   root[50];
	int    index;
	int    npart;
	double femmodel_response;
	int    flag;
	double *vertex_response   = NULL;
	double *qmu_response      = NULL;
	double *responses_pointer = NULL;

	/*retrieve npart: */
	parameters->FindParam(&npart,QmuNumberofpartitionsEnum);
	my_rank=IssmComm::GetRank();

	/*save the d_responses pointer: */
	responses_pointer=d_responses;

	//watch out, we have more d_numresponses than numresponsedescriptors, because the responses have been expanded if they were scaled. 
	//because we don't know the d_responses descriptors (the scaled ones) we can't key off them, so we will key off the responses_descriptors: */

	for(i=0;i<numresponsedescriptors;i++){

		flag=DescriptorIndex(root,&index,responses_descriptors[i]);

		if(flag==ScaledEnum){

			/*this response was scaled. pick up the response from the inputs: */
			GetVectorFromInputsx(&vertex_response,this, StringToEnumx(root),VertexPIdEnum);

			/*Now, average it onto the partition nodes: */
			AverageOntoPartitionx(&qmu_response,elements,nodes,vertices,loads,materials,parameters,vertex_response);

			/*Copy onto our dakota responses: */
			if(my_rank==0){
				/*plug response: */
				for(j=0;j<npart;j++)responses_pointer[j]=qmu_response[j];

				/*increment response_pointer :*/
				responses_pointer+=npart;
			}

			/*Free ressources:*/
			xDelete<double>(vertex_response);
			xDelete<double>(qmu_response);

		}
		else if (flag==IndexedEnum){

			/*indexed response: plug index into parameters and call response module: */
			parameters->SetParam(index,IndexEnum);

			this->Responsex(&femmodel_response,root);

			if(my_rank==0){
				/*plug response: */
				responses_pointer[0]=femmodel_response;

				/*increment response_pointer :*/
				responses_pointer++;
			}
		}
		else if (flag==NodalEnum){
			_error_("nodal response functions not supported yet!");

			/*increment response_pointer :*/
			responses_pointer++;
		}
		else if (flag==RegularEnum){

			/*perfectly normal response function: */
			this->Responsex(&femmodel_response,root);

			if(my_rank==0){
				/*plug response: */
				responses_pointer[0]=femmodel_response;

				/*increment response_pointer :*/
				responses_pointer++;
			}
		}
		else _error_("flag type " << flag << " not supported yet for response analysis");
	}

	/*Synthesize echo: {{{*/
	if(my_rank==0){
		_printf_("   responses: " << d_numresponses << ": ");
		for(i=0;i<d_numresponses-1;i++)_printf_(d_responses[i] << "|");
		_printf_(d_responses[d_numresponses-1]);
		_printf_("\n");
	}
	/*}}}*/

}
/*}}}*/
#endif
#ifdef _HAVE_GIAIVINS_
void FemModel::Deflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt, IssmDouble* x, IssmDouble* y){ /*{{{*/

	/*Go through elements, and add contribution from each element to the deflection vector wg:*/
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->GiaDeflection(wg,dwgdt, x,y);
	}
}
/*}}}*/
#endif
#ifdef _HAVE_ESA_
void FemModel::EsaGeodetic2D(Vector<IssmDouble>* pUp, Vector<IssmDouble>* pNorth, Vector<IssmDouble>* pEast, IssmDouble* xx, IssmDouble* yy){/*{{{*/

	int         ns,nsmax;
	
	/*Go through elements, and add contribution from each element to the deflection vector wg:*/
	ns = elements->Size();
	
	/*Figure out max of ns: */
	ISSM_MPI_Reduce(&ns,&nsmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(&nsmax,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Call the esa geodetic core: */
	for(int i=0;i<nsmax;i++){
		if(i<ns){
			Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
			element->EsaGeodetic2D(pUp,pNorth,pEast,xx,yy);
		}
		if(i%100==0){
			pUp->Assemble();
			pNorth->Assemble();
			pEast->Assemble();
		}
	}
	
	/*One last time: */
	pUp->Assemble();
	pNorth->Assemble();
	pEast->Assemble();

	/*Free ressources:*/
	xDelete<IssmDouble>(xx);
	xDelete<IssmDouble>(yy);
}
/*}}}*/
void FemModel::EsaGeodetic3D(Vector<IssmDouble>* pUp, Vector<IssmDouble>* pNorth, Vector<IssmDouble>* pEast, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius, IssmDouble* xx, IssmDouble* yy, IssmDouble* zz){/*{{{*/

	IssmDouble  eartharea=0;
	IssmDouble  eartharea_cpu=0;

	int         ns,nsmax;
	
	/*Go through elements, and add contribution from each element to the deflection vector wg:*/
	ns = elements->Size();
	
	/*First, figure out the surface area of Earth: */ 
	for(int i=0;i<ns;i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		eartharea_cpu += element->GetAreaSpherical();
	}
	ISSM_MPI_Reduce (&eartharea_cpu,&eartharea,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&eartharea,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Figure out max of ns: */
	ISSM_MPI_Reduce(&ns,&nsmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(&nsmax,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Call the esa geodetic core: */
	for(int i=0;i<nsmax;i++){
		if(i<ns){
			Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
			element->EsaGeodetic3D(pUp,pNorth,pEast,latitude,longitude,radius,xx,yy,zz,eartharea);
		}
		if(i%100==0){
			pUp->Assemble();
			pNorth->Assemble();
			pEast->Assemble();
		}
	}
	
	/*One last time: */
	pUp->Assemble();
	pNorth->Assemble();
	pEast->Assemble();

	/*Free ressources:*/
	xDelete<IssmDouble>(latitude);
	xDelete<IssmDouble>(longitude);
	xDelete<IssmDouble>(radius);
	xDelete<IssmDouble>(xx);
	xDelete<IssmDouble>(yy);
	xDelete<IssmDouble>(zz);
}
/*}}}*/
#endif
#ifdef _HAVE_SEALEVELRISE_
void FemModel::SealevelriseEustatic(Vector<IssmDouble>* pSgi, IssmDouble* peustatic, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius) { /*{{{*/

	/*serialized vectors:*/
	IssmDouble  eustatic       = 0.;
	IssmDouble  eustatic_cpu   = 0.;
	IssmDouble  eustatic_cpu_e = 0.;
	IssmDouble  oceanarea      = 0.;
	IssmDouble  oceanarea_cpu  = 0.;
	IssmDouble  eartharea      = 0.;
	IssmDouble  eartharea_cpu  = 0.;
	int         ns,nsmax;
	
	/*Go through elements, and add contribution from each element to the deflection vector wg:*/
	ns = elements->Size();

	/*First, figure out the area of the ocean, which is needed to compute the eustatic component: */
	for(int i=0;i<ns;i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		oceanarea_cpu += element->OceanArea();
		eartharea_cpu += element->GetAreaSpherical();
	}
	ISSM_MPI_Reduce (&oceanarea_cpu,&oceanarea,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&oceanarea,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	_assert_(oceanarea>0.);

	ISSM_MPI_Reduce (&eartharea_cpu,&eartharea,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&eartharea,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Figure out max of ns: */
	ISSM_MPI_Reduce(&ns,&nsmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(&nsmax,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Call the sea level rise core: */
	for(int i=0;i<nsmax;i++){
		if(i<ns){
		
			if(VerboseConvergence())if(i%100==0)_printf0_("\r" << "      convolution progress: " << (double)i/(double)ns*100 << "%  ");
		
			Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
			element->SealevelriseEustatic(pSgi,&eustatic_cpu_e,latitude,longitude,radius,oceanarea,eartharea);
			eustatic_cpu+=eustatic_cpu_e;
		}
		if(i%100==0)pSgi->Assemble();
	}
	if(VerboseConvergence())_printf0_("\n");
		
	/*One last time: */
	pSgi->Assemble();

	/*Sum all eustatic components from all cpus:*/
	ISSM_MPI_Reduce (&eustatic_cpu,&eustatic,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&eustatic,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	_assert_(!xIsNan<IssmDouble>(eustatic));

	/*Assign output pointers:*/
	*peustatic=eustatic;

}
/*}}}*/
void FemModel::SealevelriseNonEustatic(Vector<IssmDouble>* pSgo, Vector<IssmDouble>* pSg_old, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius, bool verboseconvolution){/*{{{*/

	/*serialized vectors:*/
	IssmDouble* Sg_old=NULL;
	
	IssmDouble  eartharea=0;
	IssmDouble  eartharea_cpu=0;

	int         ns,nsmax;
	
	/*Serialize vectors from previous iteration:*/
	Sg_old=pSg_old->ToMPISerial();

	/*Go through elements, and add contribution from each element to the deflection vector wg:*/
	ns = elements->Size();

	/*First, figure out the area of the ocean, which is needed to compute the eustatic component: */
	for(int i=0;i<ns;i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		eartharea_cpu += element->GetAreaSpherical();
	}
	
	ISSM_MPI_Reduce (&eartharea_cpu,&eartharea,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&eartharea,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Figure out max of ns: */
	ISSM_MPI_Reduce(&ns,&nsmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(&nsmax,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Call the sea level rise core: */
	for(int i=0;i<nsmax;i++){
		if(i<ns){
			if(verboseconvolution)if(VerboseConvergence())if(i%100==0)_printf_("\r" << "      convolution progress: " << (double)i/(double)ns*100 << "%   ");
			Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
			element->SealevelriseNonEustatic(pSgo,Sg_old,latitude,longitude,radius,eartharea);
		}
		if(i%100==0)pSgo->Assemble();
	}
	if(verboseconvolution)if(VerboseConvergence())_printf_("\n");
	
	/*Free ressources:*/
	xDelete<IssmDouble>(Sg_old);
}
/*}}}*/
void FemModel::SealevelriseRotationalFeedback(Vector<IssmDouble>* pSgo_rot, Vector<IssmDouble>* pSg_old, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius){/*{{{*/

	/*serialized vectors:*/
	IssmDouble* Sg_old=NULL;
	IssmDouble  eartharea=0;
	IssmDouble  eartharea_cpu=0;
	IssmDouble	tide_love_h, tide_love_k, fluid_love, moi_e, moi_p, omega, g;
	IssmDouble	load_love_k2 = -0.30922675; //degree 2 load Love number 
	IssmDouble	m1, m2, m3; 
	IssmDouble	lati, longi, radi, value; 

	/*Serialize vectors from previous iteration:*/
	Sg_old=pSg_old->ToMPISerial();

	/*First, figure out the area of the ocean, which is needed to compute the eustatic component: */
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		eartharea_cpu += element->GetAreaSpherical();
	}
	ISSM_MPI_Reduce (&eartharea_cpu,&eartharea,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&eartharea,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	IssmDouble moi_list[3]={0,0,0}; 
	IssmDouble moi_list_cpu[3]={0,0,0}; 
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->SealevelriseMomentOfInertia(&moi_list[0],Sg_old,eartharea);
		moi_list_cpu[0] += moi_list[0]; 
		moi_list_cpu[1] += moi_list[1]; 
		moi_list_cpu[2] += moi_list[2]; 
	}
	ISSM_MPI_Reduce (&moi_list_cpu[0],&moi_list[0],1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&moi_list[0],1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	// 	
	ISSM_MPI_Reduce (&moi_list_cpu[1],&moi_list[1],1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&moi_list[1],1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	// 	
	ISSM_MPI_Reduce (&moi_list_cpu[2],&moi_list[2],1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&moi_list[2],1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	
	/*pull out some useful parameters: */
	parameters->FindParam(&tide_love_h,SealevelriseTidalLoveHEnum);
	parameters->FindParam(&tide_love_k,SealevelriseTidalLoveKEnum);
	parameters->FindParam(&fluid_love,SealevelriseFluidLoveEnum);
	parameters->FindParam(&moi_e,SealevelriseEquatorialMoiEnum);
	parameters->FindParam(&moi_p,SealevelrisePolarMoiEnum);
	parameters->FindParam(&omega,SealevelriseAngularVelocityEnum);

	/*compute perturbation terms for angular velocity vector: */
	m1 = 1/(1-tide_love_k/fluid_love) * (1+load_love_k2)/(moi_p-moi_e) * moi_list[0]; 
	m2 = 1/(1-tide_love_k/fluid_love) * (1+load_love_k2)/(moi_p-moi_e) * moi_list[1]; 
	m3 = -(1+load_love_k2)/moi_p * moi_list[2];	// term associated with fluid number (3-order-of-magnitude smaller) is negelected  

	/* Green's function (1+k_2-h_2/g): checked against Glenn Milne's thesis Chapter 3 (eqs: 3.3-4, 3.10-11)
	 * Perturbation terms for angular velocity vector (m1, m2, m3): checked against Mitrovica (2005 Appendix) & Jensen et al (2013 Appendix A3) 
	 * Sea level rotational feedback: checked against GMD eqs 8-9 (only first order terms, i.e., degree 2 order 0 & 1 considered) 
	 * all DONE in Geographic coordinates: theta \in [-90,90], lambda \in [-180 180] 
	 */
	for(int i=0;i<vertices->Size();i++){
		int sid;
		//Vertex* vertex=(Vertex*)vertices->GetObjectByOffset(i);
		Vertex* vertex=xDynamicCast<Vertex*>(vertices->GetObjectByOffset(i));
		sid=vertex->Sid();

		lati=latitude[sid]/180*PI;	longi=longitude[sid]/180*PI; radi=radius[sid];

		/*only first order terms are considered now: */ 
		value=((1.0+tide_love_k-tide_love_h)/9.81)*pow(omega*radi,2.0)*
						(-m3/6.0 + 0.5*m3*cos(2.0*lati) - 0.5*sin(2.*lati)*(m1*cos(longi)+m2*sin(longi))); 
	
		pSgo_rot->SetValue(sid,value,INS_VAL); //INS_VAL ensures that you don't add several times
	}

	/*Assemble mesh velocity*/
	pSgo_rot->Assemble();
	
	/*Free ressources:*/
	xDelete<IssmDouble>(Sg_old);
	
}
/*}}}*/
void FemModel::SealevelriseGeodetic(Vector<IssmDouble>* pUp, Vector<IssmDouble>* pNorth, Vector<IssmDouble>* pEast, Vector<IssmDouble>* pSg, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius, IssmDouble* xx, IssmDouble* yy, IssmDouble* zz){/*{{{*/

	/*serialized vectors:*/
	IssmDouble* Sg=NULL;
	
	IssmDouble  eartharea=0;
	IssmDouble  eartharea_cpu=0;

	int         ns,nsmax;
	
	/*Serialize vectors from previous iteration:*/
	Sg=pSg->ToMPISerial();

	/*Go through elements, and add contribution from each element to the deflection vector wg:*/
	ns = elements->Size();
	
	/*First, figure out the area of the ocean, which is needed to compute the eustatic component: */
	for(int i=0;i<ns;i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		eartharea_cpu += element->GetAreaSpherical();
	}
	ISSM_MPI_Reduce (&eartharea_cpu,&eartharea,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&eartharea,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Figure out max of ns: */
	ISSM_MPI_Reduce(&ns,&nsmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(&nsmax,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Call the sea level rise core: */
	for(int i=0;i<nsmax;i++){
		if(i<ns){
			Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
			element->SealevelriseGeodetic(pUp,pNorth,pEast,Sg,latitude,longitude,radius,xx,yy,zz,eartharea);
		}
		if(i%100==0){
			pUp->Assemble();
			pNorth->Assemble();
			pEast->Assemble();
		}
	}
	
	/*One last time: */
	pUp->Assemble();
	pNorth->Assemble();
	pEast->Assemble();

	/*Free ressources:*/
	xDelete<IssmDouble>(Sg);
	xDelete<IssmDouble>(latitude);
	xDelete<IssmDouble>(longitude);
	xDelete<IssmDouble>(radius);
	xDelete<IssmDouble>(xx);
	xDelete<IssmDouble>(yy);
	xDelete<IssmDouble>(zz);
}
/*}}}*/
IssmDouble FemModel::SealevelriseOceanAverage(Vector<IssmDouble>* Sg) { /*{{{*/

	IssmDouble* Sg_serial=NULL;
	IssmDouble  oceanvalue,oceanvalue_cpu;
	IssmDouble  oceanarea,oceanarea_cpu;

	/*Serialize vectors from previous iteration:*/
	Sg_serial=Sg->ToMPISerial();

	/*Initialize:*/
	oceanvalue_cpu=0;
	oceanarea_cpu=0;

	/*Go through elements, and add contribution from each element and divide by overall ocean area:*/
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		oceanarea_cpu += element->OceanArea();
		oceanvalue_cpu += element->OceanAverage(Sg_serial);
	}
	ISSM_MPI_Reduce (&oceanarea_cpu,&oceanarea,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&oceanarea,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	
	ISSM_MPI_Reduce (&oceanvalue_cpu,&oceanvalue,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&oceanvalue,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Free ressources:*/
	xDelete<IssmDouble>(Sg_serial);
	
	return oceanvalue/oceanarea;
}
/*}}}*/
#endif
void FemModel::HydrologyEPLupdateDomainx(IssmDouble* pEplcount){ /*{{{*/

	Vector<IssmDouble>* mask							= NULL;
	Vector<IssmDouble>* recurence  				= NULL;
	Vector<IssmDouble>* active						= NULL;
	IssmDouble*         serial_mask				= NULL;
	IssmDouble*         serial_rec  			= NULL;
	IssmDouble*         serial_active			= NULL;
	IssmDouble*         old_active        = NULL;
	int*                eplzigzag_counter =	NULL;
	int                 eplflip_lock;
	
	HydrologyDCEfficientAnalysis* effanalysis =  new HydrologyDCEfficientAnalysis();
	HydrologyDCInefficientAnalysis* inefanalysis =  new HydrologyDCInefficientAnalysis();

	/*Step 1: update mask, the mask might be extended by residual and/or using downstream sediment head*/
	mask=new Vector<IssmDouble>(this->nodes->NumberOfNodes(HydrologyDCEfficientAnalysisEnum));
	recurence=new Vector<IssmDouble>(this->nodes->NumberOfNodes(HydrologyDCEfficientAnalysisEnum));
	this->parameters->FindParam(&eplzigzag_counter,NULL,EplZigZagCounterEnum); 
	this->parameters->FindParam(&eplflip_lock,HydrologydcEplflipLockEnum); 
	GetVectorFromInputsx(&old_active,this,HydrologydcMaskEplactiveNodeEnum,NodeSIdEnum);
	
	for (int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		effanalysis->HydrologyEPLGetMask(mask,recurence,eplzigzag_counter,element);
	}
	/*check for changes and increment zigzag counter, change the mask if necessary*/
	recurence->Assemble();
	serial_rec=recurence->ToMPISerial();
	for (int i=0;i<nodes->Size();i++){
		Node* node=xDynamicCast<Node*>(nodes->GetObjectByOffset(i));
		if(serial_rec[node->Sid()]==1.)eplzigzag_counter[node->Lid()] ++;
		if(eplzigzag_counter[node->Lid()]>eplflip_lock & eplflip_lock!=0){
			mask->SetValue(node->Sid(),old_active[node->Sid()],INS_VAL);
		}
	}
	this->parameters->SetParam(eplzigzag_counter,this->nodes->Size(),EplZigZagCounterEnum);
	/*Assemble and serialize*/
	mask->Assemble();
	serial_mask=mask->ToMPISerial();	
	
	xDelete<int>(eplzigzag_counter);
	xDelete<IssmDouble>(serial_rec);
	xDelete<IssmDouble>(old_active);
	delete mask;
	delete recurence;

	/*Update Mask*/
	InputUpdateFromVectorx(this,serial_mask,HydrologydcMaskEplactiveNodeEnum,NodeSIdEnum);
	xDelete<IssmDouble>(serial_mask);
	inefanalysis->ElementizeEplMask(this);
	/*Step 2: update node activity. If one element is connected to mask=1, all nodes are active*/
	active=new Vector<IssmDouble>(nodes->NumberOfNodes(HydrologyDCEfficientAnalysisEnum));
	for (int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		effanalysis->HydrologyEPLGetActive(active,element);
	}

	/*Assemble and serialize*/
	active->Assemble();
	serial_active=active->ToMPISerial();
	delete active;

	/*Update node activation accordingly*/
	int counter =0;
	for (int i=0;i<nodes->Size();i++){
		Node* node=xDynamicCast<Node*>(nodes->GetObjectByOffset(i));
		if(node->InAnalysis(HydrologyDCEfficientAnalysisEnum)){
			if(serial_active[node->Sid()]==1.){
				node->Activate();
				if(!node->IsClone()) counter++;
			}
			else{
				node->Deactivate();
			}
		}
	}
	xDelete<IssmDouble>(serial_active);
	delete effanalysis;
	delete inefanalysis;
	int sum_counter;
	ISSM_MPI_Reduce(&counter,&sum_counter,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&sum_counter,1,ISSM_MPI_INT,0,IssmComm::GetComm());                
	counter=sum_counter;
	*pEplcount = counter;
	if(VerboseSolution()) _printf0_("   Number of active nodes in EPL layer: "<< counter <<"\n");

	/*Update dof indexings*/
	this->UpdateConstraintsx();

}
/*}}}*/
void FemModel::UpdateConstraintsL2ProjectionEPLx(IssmDouble* pL2count){ /*{{{*/

	Vector<IssmDouble>* active        = NULL;
	IssmDouble*         serial_active = NULL;
	HydrologyDCEfficientAnalysis* effanalysis = new HydrologyDCEfficientAnalysis();

	/*update node activity. If one element is connected to mask=1, all nodes are active*/
	active=new Vector<IssmDouble>(nodes->NumberOfNodes(HydrologyDCEfficientAnalysisEnum));
	for (int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		effanalysis->HydrologyEPLGetActive(active,element);
	}

	/*Assemble and serialize*/
	active->Assemble();
	serial_active=active->ToMPISerial();
	delete active;
	delete effanalysis;

	/*Update node activation accordingly*/
	int counter =0;
	for (int i=0;i<nodes->Size();i++){
		Node* node=xDynamicCast<Node*>(nodes->GetObjectByOffset(i));
		if(node->InAnalysis(L2ProjectionEPLAnalysisEnum)){
			if(serial_active[node->Sid()]==1.){
				node->Activate();
				if(!node->IsClone()) counter++;
			}
			else{
				node->Deactivate();
			}
		}
	}
	xDelete<IssmDouble>(serial_active);
	int sum_counter;
	ISSM_MPI_Reduce(&counter,&sum_counter,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&sum_counter,1,ISSM_MPI_INT,0,IssmComm::GetComm());                
	counter=sum_counter;
	*pL2count = counter;
	if(VerboseSolution()) _printf0_("   Number of active nodes L2 Projection: "<< counter <<"\n");
}
/*}}}*/
#ifdef _HAVE_JAVASCRIPT_ 
FemModel::FemModel(IssmDouble* buffer, int buffersize, char* toolkits, char* solution, char* modelname,ISSM_MPI_Comm incomm, bool trace){ /*{{{*/
	/*configuration: */
	int  solution_type;
	int  ierr;

	/*First things first, store the communicator, and set it as a global variable: */
	IssmComm::SetComm(incomm);

	/*Start profiler: */
	this->profiler=new Profiler();
	profiler->Tag(START);

	/*From command line arguments, retrieve different filenames needed to create the FemModel: */
	solution_type=StringToEnumx(solution);
	
	/*Create femmodel from input files: */
	profiler->Tag(STARTINIT);
	this->InitFromBuffers((char*)buffer,buffersize,toolkits, solution_type,trace,NULL);
	profiler->Tag(FINISHINIT);
	
	/*Save communicator in the parameters dataset: */
	this->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(incomm,FemModelCommEnum));

}
/*}}}*/
void FemModel::CleanUpJs(char** poutput, size_t* psize){/*{{{*/

	/*Intermediary*/
	FILE *output_fid;
	GenericParam<char**>* outputbufferparam=NULL;
	GenericParam<size_t*>* outputbuffersizeparam=NULL;
	char** poutputbuffer;
	size_t* poutputbuffersize;

	
	/*Before we delete the profiler, report statistics for this run: */
	profiler->Tag(FINISH);  //final tagging
	_printf0_("\n");
	_printf0_("   "<<setw(40)<<left<<"FemModel initialization elapsed time:"<<profiler->DeltaTime(STARTINIT,FINISHINIT) << "\n");
	_printf0_("   "<<setw(40)<<left<<"Core solution elapsed time:"<<profiler->DeltaTime(STARTCORE,FINISHCORE) << "\n");
	_printf0_("\n");
	_printf0_("   Total elapsed time: "
				<<profiler->DeltaTimeModHour(START,FINISH)<<" hrs "
				<<profiler->DeltaTimeModMin(START,FINISH)<<" min "
				<<profiler->DeltaTimeModSec(START,FINISH)<<" sec"
				);
	_printf0_("\n");
	
	/*Before we close the output file, recover the buffer and size:*/
	outputbufferparam = xDynamicCast<GenericParam<char**>*>(this->parameters->FindParamObject(OutputBufferPointerEnum));
	poutputbuffer=outputbufferparam->GetParameterValue();
	outputbuffersizeparam = xDynamicCast<GenericParam<size_t*>*>(this->parameters->FindParamObject(OutputBufferSizePointerEnum));
	poutputbuffersize=outputbuffersizeparam->GetParameterValue();

	/*Assign output values: */
	*poutput=*poutputbuffer;
	*psize=*poutputbuffersize;
}
/*}}}*/
void FemModel::InitFromBuffers(char* buffer, int buffersize, char* toolkits, int in_solution_type, bool trace, IssmPDouble* X){/*{{{*/

	/*intermediary*/
	FILE       *IOMODEL = NULL;
	FILE       *toolkitsoptionsfid = NULL;
	FILE       *output_fid = NULL;
	int         my_rank;
	size_t      outputsize;
	char       *outputbuffer;
	const char *rootpath = "";   //needed for Dakota runs only, which we won't do here.

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	/*Open input file descriptor on cpu 0: */
	if(my_rank==0) IOMODEL = fmemopen((void*)buffer, buffersize, "rb");

	/*Open toolkits file descriptor: */
	toolkitsoptionsfid=fmemopen((void*)toolkits, strlen(toolkits)+1, "r");

	/*Now, go create FemModel:*/
	this->InitFromFids((char*)rootpath,IOMODEL,toolkitsoptionsfid,in_solution_type,trace,X);

	/*Close input file and toolkits file descriptors: */
	if(my_rank==0) fclose(IOMODEL);
	fclose(toolkitsoptionsfid);

	/*Open output file once for all and add output file descriptor to parameters*/
	output_fid=open_memstream(&outputbuffer,&outputsize); 
	if(output_fid==NULL)_error_("could not initialize output stream");
	this->parameters->SetParam(output_fid,OutputFilePointerEnum);
	this->parameters->AddObject(new GenericParam<char**>(&outputbuffer,OutputBufferPointerEnum));
	this->parameters->AddObject(new GenericParam<size_t*>(&outputsize,OutputBufferSizePointerEnum));

}/*}}}*/
#endif

#ifdef _HAVE_NEOPZ_
void FemModel::InitializeAdaptiveRefinement(void){/*{{{*/
	
	/*Define variables*/
	int my_rank						= IssmComm::GetRank();
	this->amr						= NULL;//initialize amr as NULL
	int numberofvertices			= this->vertices->NumberOfVertices();
	int numberofelements			= this->elements->NumberOfElements();
	int numberofsegments			= 0; //used on matlab
	IssmDouble* x					= NULL;
	IssmDouble* y					= NULL;
	IssmDouble* z					= NULL;
	int* elements					= NULL;
	int elementswidth				= this->GetElementsWidth(); //just tria elements in this version. Itapopo:
	int levelmax					= 0;
	IssmDouble regionlevel1		= 0.;
	IssmDouble regionlevelmax	= 0.;

	/*Get vertices coordinates of the coarse mesh (father mesh)*/
	/*elements comes in Matlab indexing*/
	this->GetMesh(this->vertices,this->elements,&x,&y,&z,&elements);
	
	/*Get amr parameters*/
	this->parameters->FindParam(&levelmax,AmrLevelMaxEnum);
	this->parameters->FindParam(&regionlevel1,AmrRegionLevel1Enum);
	this->parameters->FindParam(&regionlevelmax,AmrRegionLevelMaxEnum);

	/*Create initial mesh (coarse mesh) in neopz data structure*/ 
	/*Just CPU #0 should keep AMR object*/
   this->SetRefPatterns();
	if(my_rank==0){ 
	   bool ismisomip	= true;
		if(ismisomip){//itapopo
			TPZFileStream fstr;
			std::stringstream ss;

			ss							<< levelmax;
			std::string AMRfile  = "/home/santos/Misomip2/L" + ss.str() + "_tsai/amr.txt";
			fstr.OpenRead(AMRfile.c_str());
			
			TPZSaveable *sv		= TPZSaveable::Restore(fstr,0);
			this->amr				= dynamic_cast<AdaptiveMeshRefinement*>(sv);
		}
		else{
			this->amr = new AdaptiveMeshRefinement();
			//this->amr->SetLevelMax(levelmax); //Set max level of refinement
			//this->amr->SetRegions(regionlevel1,regionlevelmax);
			this->amr->CreateInitialMesh(numberofvertices,numberofelements,numberofsegments,elementswidth,x,y,z,elements,NULL);
		}
		this->amr->SetLevelMax(levelmax); //Set max level of refinement
		this->amr->SetRegions(regionlevel1,regionlevelmax);
	}

	/*Free the vectors*/
	xDelete<IssmDouble>(x);
	xDelete<IssmDouble>(y);
	xDelete<IssmDouble>(z);
	xDelete<int>(elements);

}
/*}}}*/
void FemModel::SetRefPatterns(){/*{{{*/

   /*Initialize the global variable of refinement patterns*/
   gRefDBase.InitializeUniformRefPattern(EOned);
   gRefDBase.InitializeUniformRefPattern(ETriangle);

    //gRefDBase.InitializeRefPatterns();
   /*Insert specifics patterns to ISSM core*/
   std::string filepath  = REFPATTERNDIR;
   std::string filename1 = filepath + "/2D_Triang_Rib_3.rpt";
   std::string filename2 = filepath + "/2D_Triang_Rib_4.rpt";
   std::string filename3 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_4.rpt";
   std::string filename4 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_4_permuted.rpt";
   std::string filename5 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_5.rpt";
   std::string filename6 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_5_permuted.rpt";
   std::string filename7 = filepath + "/2D_Triang_Rib_5.rpt";

   TPZAutoPointer<TPZRefPattern> refpat1 = new TPZRefPattern(filename1);
   TPZAutoPointer<TPZRefPattern> refpat2 = new TPZRefPattern(filename2);
   TPZAutoPointer<TPZRefPattern> refpat3 = new TPZRefPattern(filename3);
   TPZAutoPointer<TPZRefPattern> refpat4 = new TPZRefPattern(filename4);
   TPZAutoPointer<TPZRefPattern> refpat5 = new TPZRefPattern(filename5);
   TPZAutoPointer<TPZRefPattern> refpat6 = new TPZRefPattern(filename6);
   TPZAutoPointer<TPZRefPattern> refpat7 = new TPZRefPattern(filename7);

   if(!gRefDBase.FindRefPattern(refpat1)) gRefDBase.InsertRefPattern(refpat1);
   if(!gRefDBase.FindRefPattern(refpat2)) gRefDBase.InsertRefPattern(refpat2);
   if(!gRefDBase.FindRefPattern(refpat3)) gRefDBase.InsertRefPattern(refpat3);
   if(!gRefDBase.FindRefPattern(refpat4)) gRefDBase.InsertRefPattern(refpat4);
   if(!gRefDBase.FindRefPattern(refpat5)) gRefDBase.InsertRefPattern(refpat5);
   if(!gRefDBase.FindRefPattern(refpat6)) gRefDBase.InsertRefPattern(refpat6);
   if(!gRefDBase.FindRefPattern(refpat7)) gRefDBase.InsertRefPattern(refpat7);
}
/*}}}*/
void FemModel::ReMesh(void){/*{{{*/

	/*Variables*/
	IssmDouble *newx			= NULL;
	IssmDouble *newy			= NULL;
	IssmDouble *newz			= NULL;
	int *newelementslist		= NULL;
	int newnumberofvertices	= -1;
	int newnumberofelements = -1;
	bool* my_elements			= NULL; 
	int* my_vertices			= NULL;
	int elementswidth			= this->GetElementsWidth();//just tria elements in this version

	/*Execute refinement and get the new mesh*/
	/*newelementslist come in Matlab indexing*/
	this->ExecuteRefinement(newnumberofvertices,newnumberofelements,&newx,&newy,&newz,&newelementslist);

	/*Partitioning the new mesh. Maybe ElementsAndVerticesPartitioning.cpp could be modified to set this without iomodel.*/
	this->ElementsAndVerticesPartitioning(newnumberofvertices,newnumberofelements,elementswidth,newelementslist,&my_elements,&my_vertices);

	if(this->loads->Size()!=0) _error_("not supported yet");

	/*Create vertices*/
	Vertices* new_vertices=new Vertices();
	this->CreateVertices(newnumberofvertices,newnumberofelements,elementswidth,newelementslist,my_vertices,newx,newy,newz,new_vertices);
 
	/*Creating elements*/
	/*Just Tria in this version*/
	Elements* new_elements=new Elements();
	this->CreateElements(newnumberofelements,elementswidth,newelementslist,my_elements,new_elements);

	/*Creating materials*/
	Materials* new_materials=new Materials();
	this->CreateMaterials(newnumberofelements,my_elements,new_materials);
	
	/*Creating nodes and constraints*/
	/*Just SSA (2D) and P1 in this version*/
	Nodes* new_nodes=new Nodes();
	Constraints* new_constraints=new Constraints();

	int nodecounter=0;
	int constraintcounter=0;
	for(int i=0;i<this->nummodels;i++){//create nodes for each analysis in analysis_type_list
	
		int analysis_enum = this->analysis_type_list[i];
		
		/*As the domain is 2D, it is not necessary to create nodes for this analysis*/
		/*itapopo must verify if domain is not 3D. Only 2D in this version!*/
		if(analysis_enum==StressbalanceVerticalAnalysisEnum) continue;	    
		
		this->CreateNodes(newnumberofvertices,my_vertices,nodecounter,analysis_enum,new_nodes);
		if(analysis_enum==StressbalanceAnalysisEnum) this->CreateConstraints(newnumberofvertices,newnumberofelements,nodecounter,constraintcounter,newx,newy,my_vertices,new_constraints);
		this->UpdateElements(newnumberofelements,newelementslist,my_elements,nodecounter,i,new_elements);

		if(new_nodes->Size()) nodecounter=new_nodes->MaximumId();
		constraintcounter = new_constraints->NumberOfConstraints();
		/*Make sure nodecounter is at least 0 (if no node exists, maxid will be -1*/
		_assert_(nodecounter>=0);
	}

	new_elements->Presort();
	new_nodes->Presort();
	new_vertices->Presort();
	this->loads->Presort();
	new_materials->Presort();
	new_constraints->Presort();

	/*reset hooks for elements, loads and nodes: */
	new_elements->ResetHooks();
	this->loads->ResetHooks();
	new_materials->ResetHooks();

	/*do the post-processing of the datasets to get an FemModel that can actually run analyses: */
	int analysis_type;
	for(int i=0;i<this->nummodels;i++){
		analysis_type=this->analysis_type_list[i];
		//SetCurrentConfiguration(analysis_type);

		this->analysis_counter=i;	
		/*Now, plug analysis_counter and analysis_type inside the parameters: */
		this->parameters->SetParam(this->analysis_counter,AnalysisCounterEnum);
		this->parameters->SetParam(analysis_type,AnalysisTypeEnum);
		this->parameters->SetParam(analysis_type,ConfigurationTypeEnum);

		/*configure elements, loads and nodes, for this new analysis: */
		new_elements->SetCurrentConfiguration(new_elements,this->loads,new_nodes,new_vertices,new_materials,this->parameters);
		this->loads->SetCurrentConfiguration(new_elements,this->loads,new_nodes,new_vertices,new_materials,this->parameters);

		/*take care of toolkits options, that depend on this analysis type (present only after model processor)*/
		if(this->parameters->Exist(ToolkitsOptionsStringsEnum)){
			ToolkitsOptionsFromAnalysis(this->parameters,analysis_type);
			if(VerboseSolver()) _printf0_("      toolkits Options set for analysis type: " << EnumToStringx(analysis_type) << "\n");
		}
		
		ConfigureObjectsx(new_elements,this->loads,new_nodes,new_vertices,new_materials,this->parameters);
		if(i==0){ 
			VerticesDofx(new_vertices,this->parameters); //only call once, we only have one set of vertices
		}
		SpcNodesx(new_nodes,new_constraints,this->parameters,analysis_type);
		NodesDofx(new_nodes,this->parameters,analysis_type);
	}

	/*Finally: interpolate all inputs and insert them into the new elements.*/
	this->InterpolateInputs(new_vertices,new_elements);

	/*Delete old structure and set new pointers*/
	delete this->vertices;		this->vertices		= new_vertices;
	delete this->elements;		this->elements		= new_elements;
	delete this->nodes;			this->nodes			= new_nodes;
	delete this->constraints;	this->constraints	= new_constraints;
	delete this->materials;		this->materials	= new_materials;
	
	GetMaskOfIceVerticesLSMx(this);

	/*Insert MISMIP+ bed topography*/
	if(true) this->BedrockFromMismipPlus();
	
	/*Adjust base, thickness and mask grounded ice leve set*/
	if(true) this->AdjustBaseThicknessAndMask();

	/*Reset current configuration: */
	analysis_type=this->analysis_type_list[this->analysis_counter];
	SetCurrentConfiguration(analysis_type);

	/*Cleanup*/
	xDelete<IssmDouble>(newx);
	xDelete<IssmDouble>(newy);
	xDelete<IssmDouble>(newz);
	xDelete<int>(newelementslist);
	xDelete<int>(my_vertices);
	xDelete<bool>(my_elements);

	return;

}
/*}}}*/
void FemModel::BedrockFromMismipPlus(void){/*{{{*/

	/*Insert bedrock from mismip+ setup*/
	/*This was used to Misomip project/simulations*/
	
	IssmDouble x,y,bx,by;
	int numvertices 		= this->GetElementsWidth();
	IssmDouble* xyz_list = NULL;
   IssmDouble* r        = xNew<IssmDouble>(numvertices);

	for(int el=0;el<this->elements->Size();el++){
      Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(el));
 		
		element->GetVerticesCoordinates(&xyz_list);
		for(int i=0;i<numvertices;i++){
			x = *(xyz_list+3*i+0);
			y = *(xyz_list+3*i+1);
			
			bx=-150.-728.8*std::pow(x/300000.,2)+343.91*std::pow(x/300000.,4)-50.57*std::pow(x/300000.,6);
			by=500./(1.+std::exp((-2./4000.)*(y-80000./2.-24000.)))+500./(1.+std::exp((2./4000.)*(y-80000./2.+24000.)));
			
			r[i]=std::max(bx+by,-720.);
		}	

		/*insert new bedrock*/
		element->AddInput(BedEnum,&r[0],P1Enum);
		
		/*Cleanup*/
		xDelete<IssmDouble>(xyz_list);
	}

   /*Delete*/
   xDelete<IssmDouble>(r);
	
	return;
}
/*}}}*/
void FemModel::AdjustBaseThicknessAndMask(void){/*{{{*/

	int     numvertices = this->GetElementsWidth();
   IssmDouble rho_water,rho_ice,density,base_float;
   IssmDouble* phi     = xNew<IssmDouble>(numvertices);
   IssmDouble* h       = xNew<IssmDouble>(numvertices);
   IssmDouble* s       = xNew<IssmDouble>(numvertices);
   IssmDouble* b       = xNew<IssmDouble>(numvertices);
   IssmDouble* r       = xNew<IssmDouble>(numvertices);
   IssmDouble* sl      = xNew<IssmDouble>(numvertices);

	for(int el=0;el<this->elements->Size();el++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(el));
	
		element->GetInputListOnVertices(&s[0],SurfaceEnum);
		element->GetInputListOnVertices(&r[0],BedEnum);
		element->GetInputListOnVertices(&sl[0],SealevelEnum);
		rho_water   = element->matpar->GetMaterialParameter(MaterialsRhoSeawaterEnum);
		rho_ice     = element->matpar->GetMaterialParameter(MaterialsRhoIceEnum);
		density     = rho_ice/rho_water;

		for(int i=0;i<numvertices;i++){
			/*calculate base floatation (which supports given surface*/
			base_float = rho_ice*s[i]/(rho_ice-rho_water);
			if(r[i]>base_float){
				b[i] = r[i];			
			} 
			else {
				b[i] = base_float;	
			} 

			if(abs(sl[i])>0) _error_("Sea level value not supported!");
			/*update thickness and mask grounded ice level set*/
			h[i]	  = s[i]-b[i];
			phi[i]  = h[i]+r[i]/density;	
		}

		/*Update inputs*/
		element->AddInput(MaskGroundediceLevelsetEnum,&phi[0],P1Enum);
		element->AddInput(ThicknessEnum,&h[0],P1Enum);
		element->AddInput(BaseEnum,&b[0],P1Enum);

	}
	
   /*Delete*/
   xDelete<IssmDouble>(phi);
   xDelete<IssmDouble>(h);
   xDelete<IssmDouble>(s);
   xDelete<IssmDouble>(b);
   xDelete<IssmDouble>(r);
   xDelete<IssmDouble>(sl);

	return;
}
/*}}}*/
void FemModel::WriteMeshInResults(void){/*{{{*/

	int step					= -1;
	int numberofelements = -1;
	int numberofvertices = -1;
	IssmDouble time		= -1;
	IssmDouble* x			= NULL;
	IssmDouble* y			= NULL;
	IssmDouble* z			= NULL;
	int* elementslist		= NULL;

	if(!this->elements || !this->vertices || !this->results || !this->parameters) return;
	 
	parameters->FindParam(&step,StepEnum);
	parameters->FindParam(&time,TimeEnum);
	numberofelements=this->elements->NumberOfElements();
	numberofvertices=this->vertices->NumberOfVertices();

	/*Get mesh. Elementslist comes in Matlab indexing*/
	this->GetMesh(this->vertices,this->elements,&x,&y,&z,&elementslist);

	/*Write mesh in Results*/
	this->results->AddResult(new GenericExternalResult<int*>(this->results->Size()+1,MeshElementsEnum,
																					elementslist,numberofelements,this->GetElementsWidth(),step,time));

	this->results->AddResult(new GenericExternalResult<IssmDouble*>(this->results->Size()+1,MeshXEnum,
																					x,numberofvertices,1,step,time));

	this->results->AddResult(new GenericExternalResult<IssmDouble*>(this->results->Size()+1,MeshYEnum,
																					y,numberofvertices,1,step,time));
	
	/*Cleanup*/
	xDelete<IssmDouble>(x);
	xDelete<IssmDouble>(y);
	xDelete<IssmDouble>(z);
	xDelete<int>(elementslist);

	return;
}
/*}}}*/
void FemModel::InterpolateInputs(Vertices* newfemmodel_vertices,Elements* newfemmodel_elements){/*{{{*/

	int maxinputs = MaximumNumberOfDefinitionsEnum;

	/*Figure out how many inputs we have and their respective interpolation*/
	Vector<IssmDouble>* input_interpolations=new Vector<IssmDouble>(maxinputs);
	if(this->elements->Size()){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(0));
		element->GetInputsInterpolations(input_interpolations);
	}
	input_interpolations->Assemble();
	
	/*Serialize and set output*/
	IssmDouble* input_interpolations_serial = input_interpolations->ToMPISerial();
	delete input_interpolations;

	/*Count and get enums of all inputs in old mesh*/
	int  numP0inputs;
	int  numP1inputs;
	int *P0input_enums  = NULL;
	int *P1input_enums  = NULL;
	int *P0input_interp = NULL;
	int *P1input_interp = NULL;
	for(int step=0;step<2;step++){
		if(step){
			P0input_enums  = xNew<int>(numP0inputs);
			P1input_enums  = xNew<int>(numP1inputs);
			P0input_interp = xNew<int>(numP0inputs);
			P1input_interp = xNew<int>(numP1inputs);	
		}
		numP0inputs = 0;
		numP1inputs = 0;
		for(int i=0;i<maxinputs;i++){
			int inputinterp = reCast<int>(input_interpolations_serial[i]);
			switch(inputinterp){
				case 0:
					/*Input not found, go to the next*/
					break;
				case P1Enum:
					if(step){
						P1input_enums[numP1inputs]  = i;
						P1input_interp[numP1inputs] = inputinterp;
					}
					 numP1inputs++;
					break;
				case P0Enum:
				case DoubleInputEnum:
				case IntInputEnum:
				case BoolInputEnum:
					if(step){
						P0input_enums[numP0inputs]  = i;
						P0input_interp[numP0inputs] = inputinterp;
					}
					numP0inputs++;
					break;
				default:
					_error_(EnumToStringx(inputinterp)<<" Not supported yet");
			}
		}
	}

	/*========== Deal with P0 inputs ==========*/
	int numelementsold = this->elements->NumberOfElements();
	int numelementsnew = newfemmodel_elements->NumberOfElements();
	int numverticesold = this->vertices->NumberOfVertices();
	int numverticesnew = newfemmodel_vertices->NumberOfVertices();
	IssmDouble* P0inputsold = xNew<IssmDouble>(numelementsold*numP0inputs);
	IssmDouble* P0inputsnew = NULL;
	IssmDouble* vector      = NULL;

	for(int i=0;i<numP0inputs;i++){
		GetVectorFromInputsx(&vector,this,P0input_enums[i],ElementSIdEnum);

		/*Copy vector in matrix*/
		for(int j=0;j<numelementsold;j++) P0inputsold[j*numP0inputs+i] = vector[j];
		xDelete<IssmDouble>(vector);
	}

	/*========== Deal with P1 inputs ==========*/
	IssmDouble* P1inputsold = xNew<IssmDouble>(numverticesold*numP1inputs);
	IssmDouble* P1inputsnew = NULL;
	vector      = NULL;

	for(int i=0;i<numP1inputs;i++){
		GetVectorFromInputsx(&vector,this,P1input_enums[i],VertexSIdEnum);

		/*Copy vector in matrix*/
		for(int j=0;j<numverticesold;j++) P1inputsold[j*numP1inputs+i] = vector[j];
		xDelete<IssmDouble>(vector);
	}
	
	/*Old mesh coordinates*/
	IssmDouble *Xold     = NULL;
	IssmDouble *Yold     = NULL;
	IssmDouble *Zold		= NULL;
	int        *Indexold = NULL;
	IssmDouble *Xnew     = NULL;
	IssmDouble *Ynew     = NULL;
	IssmDouble *Znew		= NULL;
	IssmDouble* XC_new   = NULL;
	IssmDouble* YC_new   = NULL;
	int        *Indexnew = NULL;
	
	/*Get the old mesh*/
	this->GetMesh(this->vertices,this->elements,&Xold,&Yold,&Zold,&Indexold);

	/*Get the new mesh*/
	this->GetMesh(newfemmodel_vertices,newfemmodel_elements,&Xnew,&Ynew,&Znew,&Indexnew);

	/*Calculate the center points xc and xy*/
	int elementswidth = this->GetElementsWidth(); //just tria in this version
	/*new mesh*/
	XC_new=xNewZeroInit<IssmDouble>(numelementsnew);
	YC_new=xNewZeroInit<IssmDouble>(numelementsnew);
	for(int i=0;i<numelementsnew;i++){
		for(int j=0;j<elementswidth;j++){
			int vid = Indexnew[i*elementswidth+j]-1;//Transform to C indexing
			XC_new[i]+=Xnew[vid];
			YC_new[i]+=Ynew[vid];
		}
		XC_new[i]=XC_new[i]/3;
		YC_new[i]=YC_new[i]/3;
	}

	/*Interplate P0 inputs in the new mesh*/
	InterpFromMeshToMesh2dx(&P0inputsnew,Indexold,Xold,Yold,numverticesold,numelementsold,
				P0inputsold,numelementsold,numP0inputs,
				XC_new,YC_new,numelementsnew,NULL);
	
	/*Interpolate P1 inputs in the new mesh*/
	InterpFromMeshToMesh2dx(&P1inputsnew,Indexold,Xold,Yold,numverticesold,numelementsold,
				P1inputsold,numverticesold,numP1inputs,
				Xnew,Ynew,numverticesnew,NULL);
	
	/*Insert P0 inputs into the new elements.*/
	vector=NULL;
	for(int i=0;i<numP0inputs;i++){
		
		/*Get P0 input vector from the interpolated matrix*/
		vector=xNew<IssmDouble>(numelementsnew);
		for(int j=0;j<numelementsnew;j++) vector[j]=P0inputsnew[j*numP0inputs+i];//vector has values for all elements (serial)

		/*Update elements from inputs: */
		for(int j=0;j<newfemmodel_elements->Size();j++){
			Element* element=xDynamicCast<Element*>(newfemmodel_elements->GetObjectByOffset(j));
			switch(P0input_interp[i]){	
				case P0Enum:
				case DoubleInputEnum:
					element->AddInput(new DoubleInput(P0input_enums[i],vector[element->sid]));//sid because newfemmodel has just a partitioning 
					break;
				case IntInputEnum: 
					element->AddInput(new IntInput(P0input_enums[i],vector[element->sid]));//sid because newfemmodel has just a partitioning
					break;
				case BoolInputEnum:
					element->AddInput(new BoolInput(P0input_enums[i],vector[element->sid]));//sid because newfemmodel has just a partitioning
					break;
				default:
					_error_(EnumToStringx(P0input_enums[i])<<" Not supported yet");
			}
		}

		xDelete<IssmDouble>(vector);
	}

	/*Insert P1 inputs into the new elements.*/
	vector=NULL;
	for(int i=0;i<numP1inputs;i++){

		/*Get P1 input vector from the interpolated matrix*/
		vector=xNew<IssmDouble>(numverticesnew);
		for(int j=0;j<numverticesnew;j++) vector[j]=P1inputsnew[j*numP1inputs+i];//vector has all vertices	(serial)

		/*Update elements from inputs: */
		//InputUpdateFromVectorx(newfemmodel,vector,P1input_enums[i],VertexSIdEnum);//VertexSId because vector is serial in SId indexing
		for(int j=0;j<newfemmodel_elements->Size();j++){
			Element* element=xDynamicCast<Element*>(newfemmodel_elements->GetObjectByOffset(j));
			element->InputUpdateFromVector(vector,P1input_enums[i],VertexSIdEnum);
		}

		xDelete<IssmDouble>(vector);
	}

	/*Cleanup*/
	xDelete<IssmDouble>(input_interpolations_serial);
	xDelete<IssmDouble>(P0inputsold);
	xDelete<IssmDouble>(P0inputsnew);
	xDelete<int>(P0input_enums);
	xDelete<int>(P0input_interp);
	xDelete<IssmDouble>(P1inputsold);
	xDelete<IssmDouble>(P1inputsnew);
	xDelete<int>(P1input_enums);
	xDelete<int>(P1input_interp);
	xDelete<IssmDouble>(Xold);
	xDelete<IssmDouble>(Yold);
	xDelete<IssmDouble>(Zold);
	xDelete<int>(Indexold);
	xDelete<IssmDouble>(Xnew);
	xDelete<IssmDouble>(Ynew);
	xDelete<IssmDouble>(Znew);
	xDelete<IssmDouble>(XC_new);
	xDelete<IssmDouble>(YC_new);
	xDelete<int>(Indexnew);

}
/*}}}*/
void FemModel::ExecuteRefinement(int &numberofvertices,int &numberofelements,IssmDouble** px,IssmDouble** py,IssmDouble** pz,int** pelementslist){/*{{{*/
	
	/*elements is in Matlab indexing*/
	
	int my_rank					 = IssmComm::GetRank();
	int numberofsegments		 = -1;
	IssmDouble* vx				 = NULL; //itapopo this is not being used
	IssmDouble* vy				 = NULL; //itapopo this is not being used
	IssmDouble* x				 = NULL;
	IssmDouble* y				 = NULL;
	IssmDouble* z				 = NULL;
	int* elementslist			 = NULL;
	int* segments				 = NULL;
	IssmDouble* masklevelset = NULL;
   const int elementswidth  = this->GetElementsWidth();//just 2D mesh, tria elements
	
	/*Solutions which will be used to refine the elements*/
	this->GetGroundediceLevelSet(&masklevelset);//itapopo verificar se já existe um método igual a esse

	if(my_rank==0){
		int type_process=1; //1: it refines father mesh. See AdaptiveMeshRefinement.h (.cpp)
		this->amr->ExecuteRefinement(type_process,vx,vy,masklevelset,
												numberofvertices,numberofelements,numberofsegments,&x,&y,&z,&elementslist,&segments);
		if(numberofvertices<=0 || numberofelements<=0 /*|| newnumberofsegments<=0*/) _error_("Error in the refinement process.");
	}
	else{
		x=xNew<IssmDouble>(numberofvertices);
		y=xNew<IssmDouble>(numberofvertices);
		z=xNew<IssmDouble>(numberofvertices);
		elementslist=xNew<int>(numberofelements*this->GetElementsWidth());
	}

	/*Send new mesh to others CPU*/
	ISSM_MPI_Bcast(&numberofvertices,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(&numberofelements,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(x,numberofvertices,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());	
	ISSM_MPI_Bcast(y,numberofvertices,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());	
	ISSM_MPI_Bcast(z,numberofvertices,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());	
	ISSM_MPI_Bcast(elementslist,numberofelements*this->GetElementsWidth(),ISSM_MPI_INT,0,IssmComm::GetComm());	

	/*Assign the pointers*/	
	(*pelementslist) = elementslist; //Matlab indexing
	(*px)				  = x;
	(*py)				  = y;
	(*pz)				  = z;

	/*Cleanup*/
	if(segments) xDelete<int>(segments);
	xDelete<IssmDouble>(masklevelset);

}
/*}}}*/
void FemModel::GetGroundediceLevelSet(IssmDouble **pmasklevelset){/*{{{*/

	int elementswidth		= this->GetElementsWidth();//just 2D mesh, tria elements
	int numberofelements = this->elements->NumberOfElements();
	int numberofvertices = this->vertices->NumberOfVertices();

	IssmDouble* elementlevelset=xNew<IssmDouble>(elementswidth);
	int* elem_vertices=xNew<int>(elementswidth);
	Vector<IssmDouble>* vmasklevelset=new Vector<IssmDouble>(numberofvertices);

	for(int i=0;i<this->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(this->elements->GetObjectByOffset(i));
		element->GetInputListOnVertices(elementlevelset,MaskGroundediceLevelsetEnum);
		element->GetVerticesSidList(elem_vertices);
		vmasklevelset->SetValue(elem_vertices[0],elementlevelset[0],INS_VAL);
      vmasklevelset->SetValue(elem_vertices[1],elementlevelset[1],INS_VAL);
      vmasklevelset->SetValue(elem_vertices[2],elementlevelset[2],INS_VAL);
	}

   /*Assemble*/
	vmasklevelset->Assemble();
	
	/*Serialize and set output*/
	(*pmasklevelset)=vmasklevelset->ToMPISerial();

	/*Cleanup*/
	xDelete<IssmDouble>(elementlevelset);
	xDelete<int>(elem_vertices);
	delete vmasklevelset;

}
/*}}}*/
void FemModel::CreateVertices(int newnumberofvertices,int newnumberofelements,int elementswidth,int* newelementslist,int* my_vertices,IssmDouble* newx,IssmDouble* newy,IssmDouble* newz,Vertices* vertices){/*{{{*/

	/*newelementslist is in Matlab indexing*/
	
	/*Creating connectivity table*/
	int* connectivity=NULL;
	connectivity=xNewZeroInit<int>(newnumberofvertices);

	for (int i=0;i<newnumberofelements;i++){
		for (int j=0;j<elementswidth;j++){
			int vertexid = newelementslist[elementswidth*i+j];
			_assert_(vertexid>0 && vertexid-1<newnumberofvertices);//Matlab indexing
			connectivity[vertexid-1]+=1;//Matlab to C indexing
		}
	}	

	/*Create vertex and insert in vertices*/
	for(int i=0;i<newnumberofvertices;i++){
		if(my_vertices[i]){
			Vertex *newvertex=new Vertex();	
			newvertex->id=i+1;
			newvertex->sid=i;
			newvertex->pid=UNDEF;
			newvertex->x=newx[i];
			newvertex->y=newy[i];
			newvertex->z=newz[i];
			newvertex->domaintype=Domain2DhorizontalEnum;
			newvertex->sigma=0.;
			newvertex->connectivity=connectivity[i];
			newvertex->clone=false;//itapopo check this
			vertices->AddObject(newvertex);	
		} 
	}

	xDelete<int>(connectivity);
}
/*}}}*/
void FemModel::CreateElements(int newnumberofelements,int elementswidth,int* newelementslist,bool* my_elements,Elements* elements){/*{{{*/

	/*newlementslist is in Matlab indexing*/

	for(int i=0;i<newnumberofelements;i++){
		if(my_elements[i]){
			/*Create element - just tria in this version*/
			Tria *newtria=new Tria();
			newtria->id=i+1;
			newtria->sid=i;
			newtria->parameters=NULL;
			newtria->inputs=new Inputs();
			newtria->nodes=NULL;
			newtria->vertices=NULL;
			newtria->material=NULL;
			newtria->matpar=NULL;
			if(this->nummodels>0){
				newtria->element_type_list=xNew<int>(this->nummodels);
				for(int j=0;j<nummodels;j++) newtria->element_type_list[j]=0;
			}
			else newtria->element_type_list=NULL;
			
			/*Element hook*/
			int matpar_id=newnumberofelements+1; //retrieve material parameter id (last pointer in femodel->materials)
			int material_id=i+1; // retrieve material_id = i+1;
			/*retrieve vertices ids*/
			int* vertex_ids=xNew<int>(elementswidth);
			for(int j=0;j<elementswidth;j++)	vertex_ids[j]=reCast<int>(newelementslist[elementswidth*i+j]);//this Hook wants Matlab indexing	
			/*Setting the hooks*/
			newtria->numanalyses =this->nummodels;
			newtria->hnodes		=new Hook*[this->nummodels];
			newtria->hvertices   =new Hook(&vertex_ids[0],elementswidth);
			newtria->hmaterial   =new Hook(&material_id,1);
			newtria->hmatpar     =new Hook(&matpar_id,1);
			newtria->hneighbors  =NULL;
			/*Initialize hnodes as NULL*/
			for(int j=0;j<this->nummodels;j++) newtria->hnodes[j]=NULL;
			/*Clean up*/
			xDelete<int>(vertex_ids);
			elements->AddObject(newtria);	
		} 
	}

}
/*}}}*/
void FemModel::CreateMaterials(int newnumberofelements,bool* my_elements,Materials* materials){/*{{{*/

	/*Just Matice in this version*/
	for(int i=0;i<newnumberofelements;i++){
		if(my_elements[i]){
			materials->AddObject(new Matice(i+1,i,MaticeEnum));	
		} 
	}
	
	/*Add new constant material property to materials, at the end: */
	Matpar *newmatpar=static_cast<Matpar*>(this->materials->GetObjectByOffset(this->materials->Size()-1)->copy());
	newmatpar->SetMid(newnumberofelements+1);
	materials->AddObject(newmatpar);//put it at the end of the materials	    

}
/*}}}*/
void FemModel::CreateNodes(int newnumberofvertices,int* my_vertices,int nodecounter,int analysis_enum,Nodes* nodes){/*{{{*/

	int lid=0;
	for(int j=0;j<newnumberofvertices;j++){
		if(my_vertices[j]){				
			
			Node* newnode=new Node();	
			
			/*id: */
			newnode->id=nodecounter+j+1;
			newnode->sid=j;
			newnode->lid=lid++;
			newnode->analysis_enum=analysis_enum;
			
			/*Initialize coord_system: Identity matrix by default*/
			for(int k=0;k<3;k++) for(int l=0;l<3;l++) newnode->coord_system[k][l]=0.0;
			for(int k=0;k<3;k++) newnode->coord_system[k][k]=1.0;
			
			/*indexing:*/
			newnode->indexingupdate=true;
			
			Analysis* analysis=EnumToAnalysis(analysis_enum);
			int *doftypes=NULL;
			int numdofs=analysis->DofsPerNode(&doftypes,Domain2DhorizontalEnum,SSAApproximationEnum);
			newnode->indexing.Init(numdofs,doftypes);
			xDelete<int>(doftypes);
			delete analysis;
			if(analysis_enum==StressbalanceAnalysisEnum)
				newnode->SetApproximation(SSAApproximationEnum);
			else
				newnode->SetApproximation(0);

			/*Stressbalance Horiz*/
			if(analysis_enum==StressbalanceAnalysisEnum){
				// itapopo this code is rarely used. 
				/*Coordinate system provided, convert to coord_system matrix*/
				//XZvectorsToCoordinateSystem(&this->coord_system[0][0],&iomodel->Data(StressbalanceReferentialEnum)[j*6]);
				//_assert_(sqrt( coord_system[0][0]*coord_system[0][0] + coord_system[1][0]*coord_system[1][0]) >1.e-4);

			}
			nodes->AddObject(newnode);
		}
	}
	return;
}
/*}}}*/
void FemModel::GetMesh(Vertices* femmodel_vertices, Elements* femmodel_elements,IssmDouble** px, IssmDouble** py, IssmDouble** pz, int** pelementslist){/*{{{*/

	if(!femmodel_vertices) _error_("GetMesh: vertices are NULL.");
	if(!femmodel_elements) _error_("GetMesh: elements are NULL.");
	
	int numberofvertices, numberofelements;
	int elementswidth = this->GetElementsWidth(); // just 2D mesh in this version (just tria elements)
	IssmDouble *x		= NULL;
	IssmDouble *y		= NULL;
	IssmDouble *z		= NULL;	
	int* elementslist = NULL;
	
	/*Get vertices coordinates*/
	VertexCoordinatesx(&x, &y, &z, femmodel_vertices,false) ;

	numberofvertices = femmodel_vertices->NumberOfVertices();
	numberofelements = femmodel_elements->NumberOfElements();
	
	/*Get element vertices*/
	int* elem_vertices=xNew<int>(elementswidth);
	Vector<IssmDouble>* vid1= new Vector<IssmDouble>(numberofelements);
	Vector<IssmDouble>* vid2= new Vector<IssmDouble>(numberofelements);
	Vector<IssmDouble>* vid3= new Vector<IssmDouble>(numberofelements);

	/*Go through elements, and for each element, get vertices*/
   for(int i=0;i<femmodel_elements->Size();i++){
    	Element* element=xDynamicCast<Element*>(femmodel_elements->GetObjectByOffset(i));
    	element->GetVerticesSidList(elem_vertices);
    	vid1->SetValue(element->sid,elem_vertices[0],INS_VAL);
    	vid2->SetValue(element->sid,elem_vertices[1],INS_VAL);
    	vid3->SetValue(element->sid,elem_vertices[2],INS_VAL);
   }
		
	/*Assemble*/
   vid1->Assemble();
   vid2->Assemble();
   vid3->Assemble();

   /*Serialize*/
	IssmDouble *id1 = vid1->ToMPISerial();
   IssmDouble *id2 = vid2->ToMPISerial();
	IssmDouble *id3 = vid3->ToMPISerial();
	
	/*Construct elements list*/
	elementslist=xNew<int>(numberofelements*elementswidth);
	if(numberofelements*elementswidth<0) _error_("numberofelements negative.");
	for(int i=0;i<numberofelements;i++){
		elementslist[elementswidth*i+0] = (int)id1[i]+1; //InterpMesh wants Matlab indexing
		elementslist[elementswidth*i+1] = (int)id2[i]+1; //InterpMesh wants Matlab indexing
		elementslist[elementswidth*i+2] = (int)id3[i]+1; //InterpMesh wants Matlab indexinf
	}
	
	/*Assign pointers*/
	*px				= x;
	*py				= y;
	*pz				= z;
	*pelementslist = elementslist; //Matlab indexing. InterMesh uses this type.

	/*Cleanup*/
	xDelete<int>(elem_vertices);
	xDelete<IssmDouble>(id1);
	xDelete<IssmDouble>(id2);
	xDelete<IssmDouble>(id3);
	delete vid1;
	delete vid2;
	delete vid3;
}
/*}}}*/
void FemModel::CreateConstraints(int newnumberofvertices,int newnumberofelements,int nodecounter,int constraintcounter,IssmDouble* newx,IssmDouble* newy,int* my_vertices,Constraints* constraints){/*{{{*/

	/*itapopo ATTENTION: JUST SPCVX AND SPCVY TO TEST!!!*/
	/*OTHERS CONSTRAINTS MUST BE IMPLEMENTED!!!*/
	
	/*Get x and y of the mesh i-1*/
	int numberofvertices			= this->vertices->NumberOfVertices();
	int numberofelements			= this->elements->NumberOfElements();
	IssmDouble *x					= NULL;
	IssmDouble *y					= NULL;
	IssmDouble *z					= NULL;
	int        *elementslist	= NULL;

	/*elementslist is in Matlab indexing*/
	this->GetMesh(this->vertices,this->elements,&x,&y,&z,&elementslist);
	
	/*Get spcvx and spcvy for mesh i-1*/
	IssmDouble *spcvx						= NULL;
	IssmDouble *spcvy						= NULL;
	IssmDouble *spcvxflag				= NULL;
	IssmDouble *spcvyflag				= NULL;
	int numberofnodes_analysistype	= this->nodes->NumberOfNodes(StressbalanceAnalysisEnum);
	Vector<IssmDouble>* vspcvx			= new Vector<IssmDouble>(numberofnodes_analysistype);
	Vector<IssmDouble>* vspcvy			= new Vector<IssmDouble>(numberofnodes_analysistype);
	Vector<IssmDouble>* vspcvxflag	= new Vector<IssmDouble>(numberofnodes_analysistype);
	Vector<IssmDouble>* vspcvyflag	= new Vector<IssmDouble>(numberofnodes_analysistype);
	
	for(int i=0;i<this->constraints->Size();i++){
		SpcStatic* spc			= xDynamicCast<SpcStatic*>(this->constraints->GetObjectByOffset(i));
		int dof					= spc->GetDof();
		int node					= spc->GetNodeId();
		IssmDouble spcvalue	= spc->GetValue(); 
		int nodeindex			= node-1;
		
		if(dof==0) {//vx
			vspcvx->SetValue(nodeindex,spcvalue,INS_VAL);
			vspcvxflag->SetValue(nodeindex,1,INS_VAL);
		} 
		if(dof==1){//vy
			vspcvy->SetValue(nodeindex,spcvalue,INS_VAL);
			vspcvyflag->SetValue(nodeindex,1,INS_VAL);
		}
	}

	/*Assemble*/
	vspcvx->Assemble();
	vspcvy->Assemble();
	vspcvxflag->Assemble();
	vspcvyflag->Assemble();

	/*Serialize*/
	spcvx		 = vspcvx->ToMPISerial();
	spcvy		 = vspcvy->ToMPISerial();
	spcvxflag = vspcvxflag->ToMPISerial();
	spcvyflag = vspcvyflag->ToMPISerial();
	
	IssmDouble *newspcvx		 = NULL;
	IssmDouble *newspcvy		 = NULL;
	IssmDouble *newspcvxflag = NULL;
	IssmDouble *newspcvyflag = NULL;
	int nods_data				 = numberofvertices;
	int nels_data				 = numberofelements;
	int M_data					 = numberofvertices;
	int N_data					 = 1;
	int N_interp				 = newnumberofvertices;

  /*Interpolate spcvx and spcvy in the new mesh*/
	InterpFromMeshToMesh2dx(&newspcvx,elementslist,x,y,nods_data,nels_data,spcvx,M_data,N_data,newx,newy,N_interp,NULL);
	InterpFromMeshToMesh2dx(&newspcvy,elementslist,x,y,nods_data,nels_data,spcvy,M_data,N_data,newx,newy,N_interp,NULL);
	InterpFromMeshToMesh2dx(&newspcvxflag,elementslist,x,y,nods_data,nels_data,spcvxflag,M_data,N_data,newx,newy,N_interp,NULL);
	InterpFromMeshToMesh2dx(&newspcvyflag,elementslist,x,y,nods_data,nels_data,spcvyflag,M_data,N_data,newx,newy,N_interp,NULL);

	int count					= 0;
	IssmDouble eps				= 1.e-8;

	/*Now, insert the interpolated constraints in the data set (constraints)*/
	for(int i=0;i<newnumberofvertices;i++){
		if(my_vertices[i]){ 
			/*spcvx*/
			if(!xIsNan<IssmDouble>(newspcvx[i]) && newspcvxflag[i]>(1-eps)){
				constraints->AddObject(new SpcStatic(constraintcounter+count+1,nodecounter+i+1,0,newspcvx[i],StressbalanceAnalysisEnum));
				//add count'th spc, on node i+1, setting dof 1 to vx.
				count++;
			}
		}
	}
	count=0;
	for(int i=0;i<newnumberofvertices;i++){
		if(my_vertices[i]){
			/*spcvy*/
			if(!xIsNan<IssmDouble>(newspcvy[i]) && newspcvyflag[i]>(1-eps) ){
				constraints->AddObject(new SpcStatic(constraintcounter+count+1,nodecounter+i+1,1,newspcvy[i],StressbalanceAnalysisEnum)); 
				//add count'th spc, on node i+1, setting dof 1 to vx.
				count++;
			}
		}
	}

	/*Cleanup*/
	xDelete<IssmDouble>(x);
	xDelete<IssmDouble>(y);
	xDelete<IssmDouble>(z);
	xDelete<int>(elementslist);
	xDelete<IssmDouble>(spcvx);
	xDelete<IssmDouble>(spcvy);
	xDelete<IssmDouble>(spcvxflag);
	xDelete<IssmDouble>(spcvyflag);
	xDelete<IssmDouble>(newspcvx);
	xDelete<IssmDouble>(newspcvy);
	xDelete<IssmDouble>(newspcvxflag);
	xDelete<IssmDouble>(newspcvyflag);
	delete vspcvx;
	delete vspcvy;	
	delete vspcvxflag;
	delete vspcvyflag;

}
/*}}}*/
void FemModel::UpdateElements(int newnumberofelements,int* newelementslist,bool* my_elements,int nodecounter,int analysis_counter,Elements* newelements){/*{{{*/

	/*newelementslist is in Matlab indexing*/

	/*Update elements, set hnode.
	This code is in all analysis */
	int elemcounter=0;
	for(int iel=0;iel<newnumberofelements;iel++){
		if(my_elements[iel]){
			Tria* tria=(Tria*)newelements->GetObjectByOffset(elemcounter);
			//element update
			tria->element_type_list[analysis_counter]=P1Enum;
			int numnodes=3;
         int* tria_node_ids=xNew<int>(numnodes);
         tria_node_ids[0]=nodecounter+newelementslist[3*iel+0]; //matlab indexing
         tria_node_ids[1]=nodecounter+newelementslist[3*iel+1]; //matlab indexing
         tria_node_ids[2]=nodecounter+newelementslist[3*iel+2]; //matlab indexing
			tria->SetHookNodes(tria_node_ids,numnodes,analysis_counter); tria->nodes=NULL;
   		xDelete<int>(tria_node_ids);
			elemcounter++;
		}
	}
	return;
}
/*}}}*/
void FemModel::ElementsAndVerticesPartitioning(int& newnumberofvertices,int& newnumberofelements,int& elementswidth,int* newelementslist,bool** pmy_elements,int** pmy_vertices){/*{{{*/

	/*newelementslist come in Matlab indexing*/

	int *epart			= NULL; //element partitioning.
	int *npart			= NULL; //node partitioning.
	int *index 			= NULL; //elements in C indexing
	int edgecut			= 1;
	int numflag			= 0;
	int etype			= 1;
	int my_rank			= IssmComm::GetRank();
	int numprocs		= IssmComm::GetSize();
	bool *my_elements = NULL;
	int *my_vertices  = NULL;
	
	_assert_(newnumberofvertices>0); 
	_assert_(newnumberofelements>0); 
	epart=xNew<int>(newnumberofelements);
	npart=xNew<int>(newnumberofvertices);
   index=xNew<int>(elementswidth*newnumberofelements);
   
	for (int i=0;i<newnumberofelements;i++){
   	for (int j=0;j<elementswidth;j++){
      	*(index+elementswidth*i+j)=(*(newelementslist+elementswidth*i+j))-1; //-1 for C indexing in Metis
      }
   }

	/*Partition using Metis:*/
	if (numprocs>1){
#ifdef _HAVE_METIS_
		METIS_PartMeshNodalPatch(&newnumberofelements,&newnumberofvertices, index, &etype, &numflag, &numprocs, &edgecut, epart, npart);
#else
		_error_("metis has not beed installed. Cannot run with more than 1 cpu");
#endif
	}
	else if (numprocs==1){
		/*METIS does not know how to deal with one cpu only!*/
		for (int i=0;i<newnumberofelements;i++) epart[i]=0;
		for (int i=0;i<newnumberofvertices;i++) npart[i]=0;
	}
	else _error_("At least one processor is required");	    

	my_vertices=xNew<int>(newnumberofvertices);
	my_elements=xNew<bool>(newnumberofelements);
	for(int i=0;i<newnumberofvertices;i++) my_vertices[i]=0;
	for(int i=0;i<newnumberofelements;i++) my_elements[i]=false;

	/*Start figuring out, out of the partition, which elements belong to this cpu: */
	for(int i=0;i<newnumberofelements;i++){
		/*!All elements have been partitioned above, only deal with elements for this cpu: */
		if(my_rank==epart[i]){ 
			my_elements[i]=true;
			/*Now that we are here, we can also start building the list of vertices belonging to this cpu partition: we use 
			 *the  element index to do this. For each element n, we know index[n][0:2] holds the indices (matlab indexing) 
			 into the vertices coordinates. If we start plugging 1 into my_vertices for each index[n][i] (i=0:2), then my_vertices 
			 will hold which vertices belong to this partition*/
			for(int j=0;j<elementswidth;j++){
				_assert_(newelementslist[elementswidth*i+j]-1<newnumberofvertices);//newelementslist is in Matlab indexing 
				my_vertices[newelementslist[elementswidth*i+j]-1]=1;//newelementslist is in Matlab indexing
			}
		}
	}

	/*Assign output pointers:*/
	*pmy_elements=my_elements;
	*pmy_vertices=my_vertices;

	/*Free ressources:*/
	xDelete<int>(epart);
	xDelete<int>(npart);	    
	xDelete<int>(index);

}
/*}}}*/
#endif
