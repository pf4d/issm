/*!\file:  issm_ocean.cpp
 * \brief: ISSM OCEAN main program. 
 */ 

#include "./issm.h"

int main(int argc,char **argv){

//	/*diverse:*/
//	int    nummodels = 2;
//	int*   commsizes=xNew<int>(nummodels);
//	int*   rankzeros=xNew<int>(nummodels);
//	char** modelnames=xNew<char*>(nummodels);
//	char** dirnames=xNew<char*>(nummodels);
//	int    iceid=0; 
//	int    oceanid=1; 
//	int    modelid;
//	int    my_rank;
//	int    count=0;
//	ISSM_MPI_Comm worldcomm;
//	ISSM_MPI_Comm modelcomm;
//	ISSM_MPI_Comm toonceancomm;
//	ISSM_MPI_Comm fromicecomms;

	/*Initialize exception trapping: */
	ExceptionTrapBegin();

//	/*Initialize environment (MPI, PETSC, MUMPS, etc ...)*/
//	worldcomm=EnvironmentInit(argc,argv);
//	
//	/*What is my rank?:*/
//	ISSM_MPI_Comm_rank(worldcomm,&my_rank);
//
//	/*First model is ice, second is ocean*/
//	for(int i=0;i<nummodels;i++){
//		char* string=NULL;
//		
//		string=xNew<char>(strlen(argv[5+3*i])+1);
//		xMemCpy<char>(string,argv[5+3*i],strlen(argv[5+3*i])+1);
//		dirnames[i]=string;
//		
//		string=xNew<char>(strlen(argv[5+3*i+1])+1);
//		xMemCpy<char>(string,argv[5+3*i+1],strlen(argv[5+3*i+1])+1);
//		modelnames[i]=string;
//
//		commsizes[i]=(int) strtol(argv[5+3*i+2], (char **)NULL, 10);
//	}
//
//	/*Figure out which model each cpu will belong to: */
//	count=0;
//	for(int i=0;i<nummodels;i++){
//		if(my_rank>=count && my_rank<(count+commsizes[i])){
//			modelid=i;
//			break;
//		}
//		count+=commsizes[i];
//	} 
//	/*Buil array of who is rank 0 of their own group:*/
//	count=0;
//	for(int i=0;i<nummodels;i++){
//		rankzeros[i]=count;
//		count+=commsizes[i];
//	}
//	/*}}}*/
//
//	/*Split world into sub-communicators for each and every model:*/
//	ISSM_MPI_Comm_split(worldcomm,modelid, my_rank, &modelcomm);
//
//	/*Build inter communicators:*/ // change to Dimitris solution
//	if(modelid==iceid){
//		ISSM_MPI_Intercomm_create( modelcomm, 0, worldcomm, rankzeros[oceanid], iceid, fromicecomms+i); //communicate from local erth comm 9rank 0) to ice comm (rank 0) using modelid tag.
//	}
//	else{
//		ISSM_MPI_Intercomm_create( modelcomm, 0, worldcomm, rankzeros[iceid], oceanid, &toearthcomm); //communicate from local ice comm (rank 0) to earth comm (rank 0) using modelid tag.
//	}

//	/*Supply specific argc and argv for each sub-communicator (corresponding to each  model specificatiions):{{{*/
//	char** arguments=xNew<char*>(4);
//	arguments[0]=xNew<char>(strlen(argv[0])+1); xMemCpy<char>(arguments[0],argv[0],strlen(argv[0])+1); //executable name
//	arguments[1]=xNew<char>(strlen(argv[1])+1); xMemCpy<char>(arguments[1],argv[1],strlen(argv[1])+1); //solution name
//	arguments[2]=xNew<char>(strlen(argv[5+3*modelid])+1); xMemCpy<char>(arguments[2],argv[5+3*modelid],strlen(argv[5+3*modelid])+1); //directory name
//	arguments[3]=xNew<char>(strlen(argv[5+3*modelid+1])+1); xMemCpy<char>(arguments[3],argv[5+3*modelid+1],strlen(argv[5+3*modelid+1])+1); //model name
//	/*}}}*/
//

	//REMOVE
	/*Initialize environment (MPI, PETSC, MUMPS, etc ...)*/
	ISSM_MPI_Comm comm_init=EnvironmentInit(argc,argv);
	/*Initialize femmodel from arguments provided command line: */
	FemModel *femmodel = new FemModel(argc,argv,comm_init);
	///*Initialize femmodel from arguments provided command line: */
	//FemModel *femmodel = new FemModel(1,arguments,modelcomm);
	
//	/*Now that the models are initialized, keep communicator information in the parameters datasets of each model: */
//	femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(worldcomm,WorldCommEnum));
//	femmodel->parameters->AddObject(new IntParam(NumModelsEnum,nummodels));
//	femmodel->parameters->AddObject(new IntParam(ModelIdEnum,oceanid));
//	femmodel->parameters->AddObject(new IntParam(EarthIdEnum,iceid));
//	if(modelid==earthid) femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm*>(fromicecomms,IcecapToEarthCommEnum));
//	else femmodel->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(toearthcomm,IcecapToEarthCommEnum));

	/*Solve: */
	femmodel->Solve();

	/*Output results: */
	OutputResultsx(femmodel);

	/*Wrap up: */
	femmodel->CleanUp();

	/*Delete Model: */
	delete femmodel;

	/*Finalize environment:*/
	EnvironmentFinalize();

	/*Finalize exception trapping: */
	ExceptionTrapEnd();

//	/*Free ressources:*/
//	xDelete<int>(commsizes);
//	for(int i=0;i<nummodels;i++){
//		char* string=NULL;
//		string=dirnames[i]; xDelete<char>(string);
//		string=modelnames[i]; xDelete<char>(string);
//	}
//	xDelete<char*>(dirnames);
//	xDelete<char*>(modelnames);

	/*Return unix success: */
	return 0; 
}
