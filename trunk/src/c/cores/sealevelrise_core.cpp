/*!\file: sealevelrise_core.cpp
 * \brief: core of the SLR solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void sealevelrise_core(FemModel* femmodel){ /*{{{*/

	Vector<IssmDouble> *Sg    = NULL;
	Vector<IssmDouble> *Sg_absolute  = NULL; 
	Vector<IssmDouble> *Sg_eustatic  = NULL; 
	Vector<IssmDouble> *U_radial  = NULL; 
	Vector<IssmDouble> *U_north   = NULL; 
	Vector<IssmDouble> *U_east    = NULL; 
	bool save_results,isslr,iscoupler;
	int configuration_type;
	int solution_type;
	int        numoutputs        = 0;
	char     **requested_outputs = NULL;
	
	/*additional parameters: */
	int  gsize;
	bool spherical=true;
	IssmDouble          *latitude   = NULL;
	IssmDouble          *longitude  = NULL;
	IssmDouble          *radius     = NULL;
	IssmDouble          *xx     = NULL;
	IssmDouble          *yy     = NULL;
	IssmDouble          *zz     = NULL;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&isslr,TransientIsslrEnum);
	femmodel->parameters->FindParam(&iscoupler,TransientIscouplerEnum);

	/*first, recover lat,long and radius vectors from vertices: */
	VertexCoordinatesx(&latitude,&longitude,&radius,femmodel->vertices,spherical); 

	/*recover x,y,z vectors from vertices: */
	VertexCoordinatesx(&xx,&yy,&zz,femmodel->vertices); 

	/*Figure out size of g-set deflection vector and allocate solution vector: */
	gsize      = femmodel->nodes->NumberOfDofs(configuration_type,GsetEnum);

	/*several cases here, depending on value of iscoupler and isslr: 
	solution_type == SealevelriseSolutionEnum)       we are running sea level rise core (no coupler)
	( !iscoupler & !isslr)       we are not interested in being here :) 
	( !iscoupler & isslr)        we are running in uncoupled mode
	( iscoupler & isslr)         we are running in coupled mode, and better be earth
	( iscoupler & !isslr)        we are running in coupled mode, and better be an ice cap
	*/

	if(solution_type==SealevelriseSolutionEnum){
		isslr=1;
		iscoupler=0;
	}

	/*early return: */
	if( !iscoupler & !isslr) return;  //we are not interested in being here :) 

	/*In what follows we assume we are all running slr, either in coupled, or uncoupled mode:*/
	if(VerboseSolution()) _printf0_("   computing sea level rise\n");

	/*set configuration: */
	if(isslr)femmodel->SetCurrentConfiguration(SealevelriseAnalysisEnum);

	/*transfer deltathickness forcing from ice caps to earth model: */
	if(iscoupler) TransferForcing(femmodel,SealevelriseDeltathicknessEnum);

	/*call sea-level rise sub cores:*/
	if(isslr){
		Sg_eustatic=sealevelrise_core_eustatic(femmodel); //generalized eustatic (Farrel and Clark, Eq 4, 1st, 3rd and 4rd terms on the RHS.

		Sg=sealevelrise_core_noneustatic(femmodel,Sg_eustatic); //ocean loading tems  (2nd and 5th terms on the RHS of Farrel and Clark)
		
		/*get results into elements:*/
		//InputUpdateFromSolutionx(femmodel,Sg);		// from Eric 
		InputUpdateFromVectorx(femmodel,Sg,SealevelEnum,VertexSIdEnum);

		/*compute other geodetic signatures, such as absolute sea level chagne, components of 3-D crustal motion: */
		/*Initialize:*/
		U_radial = new Vector<IssmDouble>(gsize);
		U_north = new Vector<IssmDouble>(gsize);
		U_east = new Vector<IssmDouble>(gsize);
		Sg_absolute = new Vector<IssmDouble>(gsize); 
		
		/*call the geodetic main modlule:*/ 
		femmodel->SealevelriseGeodetic(U_radial,U_north,U_east,Sg,latitude,longitude,radius,xx,yy,zz); 

		/*compute: absolute sea level change = relative sea level change + vertical motion*/
		Sg->Copy(Sg_absolute); Sg_absolute->AXPY(U_radial,1); 
		
		/*get results into elements:*/
		InputUpdateFromVectorx(femmodel,U_radial,SealevelUmotionEnum,VertexSIdEnum);	// radial displacement 
		InputUpdateFromVectorx(femmodel,U_north,SealevelNmotionEnum,VertexSIdEnum);	// north motion 
		InputUpdateFromVectorx(femmodel,U_east,SealevelEmotionEnum,VertexSIdEnum);		// east motion 
		InputUpdateFromVectorx(femmodel,Sg_absolute,SealevelAbsoluteEnum,VertexSIdEnum);
		
		if(save_results){
			if(VerboseSolution()) _printf0_("   saving results\n");
			femmodel->parameters->FindParam(&requested_outputs,&numoutputs,SealevelriseRequestedOutputsEnum);
			femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
		}
			
		if(solution_type==SealevelriseSolutionEnum)femmodel->RequestedDependentsx();

		/*Free ressources:*/	
		delete Sg;
		delete Sg_eustatic;
		delete U_radial;
		delete U_north;
		delete U_east;
		delete Sg_absolute;
		if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}
	}

	/*transfer sea-level back to ice caps: */
	if(iscoupler)TransferSealevel(femmodel,SealevelEnum);
} 
/*}}}*/
void TransferForcing(FemModel* femmodel,int forcingenum){ /*{{{*/

	/*forcing being transferred from models to earth: */
	IssmDouble** forcings=NULL;
	IssmDouble*  forcing=NULL; 
	Vector<IssmDouble>* forcingglobal=NULL; 
	int*         nvs=NULL;
	
	/*transition vectors:*/
	IssmDouble** transitions=NULL;
	int          ntransitions; 
	int*         transitions_m=NULL;
	int*         transitions_n=NULL;
	int          nv;
	
	/*communicators:*/
	ISSM_MPI_Comm tocomm;
	ISSM_MPI_Comm* fromcomms=NULL;
	ISSM_MPI_Status status;
	int         my_rank;
	int         modelid,earthid;
	int         nummodels;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&modelid,ModelIdEnum);
	femmodel->parameters->FindParam(&earthid,EarthIdEnum);
	femmodel->parameters->FindParam(&nummodels,NumModelsEnum);
	my_rank=IssmComm::GetRank();
	
	/*retrieve the inter communicators that will be used to send data from each ice cap to the earth: */
	if(modelid==earthid){
		GenericParam<ISSM_MPI_Comm*>* parcoms = dynamic_cast<GenericParam<ISSM_MPI_Comm*>*>(femmodel->parameters->FindParamObject(IcecapToEarthCommEnum));
		if(!parcoms)_error_("TransferForcing error message: could not find IcecapToEarthComm communicator");
		fromcomms=parcoms->GetParameterValue();
	}
	else {
		GenericParam<ISSM_MPI_Comm>* parcom = dynamic_cast<GenericParam<ISSM_MPI_Comm>*>(femmodel->parameters->FindParamObject(IcecapToEarthCommEnum));
		if(!parcom)_error_("TransferForcing error message: could not find IcecapToEarthComm communicator");
		tocomm=parcom->GetParameterValue();
	}

	/*For each icecap, retrieve the forcing vector that will be sent to the earth model: */
	if(modelid!=earthid){
		nv=femmodel->vertices->NumberOfVertices();
		GetVectorFromInputsx(&forcing,femmodel,forcingenum,VertexSIdEnum);
	}

	/*Send the forcing to the earth model:{{{*/
	if(my_rank==0){
		if(modelid==earthid){
			forcings=xNew<IssmDouble*>(nummodels-1);
			nvs=xNew<int>(nummodels-1);
			for(int i=0;i<earthid;i++){
				ISSM_MPI_Recv(nvs+i, 1, ISSM_MPI_INT, 0,i, fromcomms[i], &status);
				forcings[i]=xNew<IssmDouble>(nvs[i]);
				ISSM_MPI_Recv(forcings[i], nvs[i], ISSM_MPI_DOUBLE, 0,i, fromcomms[i], &status);
			}
			
		}
		else{
			ISSM_MPI_Send(&nv, 1, ISSM_MPI_INT, 0, modelid, tocomm);
			ISSM_MPI_Send(forcing, nv, ISSM_MPI_DOUBLE, 0, modelid, tocomm);
		}
	}
	/*}}}*/

	/*On the earth model, consolidate all the forcings into one, and update the elements dataset accordingly: {{{*/
	if(modelid==earthid){
		
		/*Out of all the delta thicknesses, build one delta thickness vector made of all the ice cap contributions. 
		 *First, build the global delta thickness vector in the earth model: */
		nv=femmodel->vertices->NumberOfVertices();
		forcingglobal= new Vector<IssmDouble>(nv);

		/*Retrieve transition vectors, used to plug from each ice cap into the global forcing:*/
		femmodel->parameters->FindParam(&transitions,&ntransitions,&transitions_m,&transitions_n,SealevelriseTransitionsEnum);

		if(ntransitions!=earthid)_error_("TransferForcing error message: number of transition vectors is not equal to the number of icecaps!");

		/*Go through all the delta thicknesses coming from each ice cap: */
		if(my_rank==0){
			for(int i=0;i<earthid;i++){

				IssmDouble* forcingfromcap= forcings[i]; //careful, this only exists on rank 0 of the earth model!
				IssmDouble* transition=transitions[i];
				int         M=transitions_m[i];

				/*build index to plug values: */
				int*        index=xNew<int>(M); for(int i=0;i<M;i++)index[i]=reCast<int>(transition[i])-1; //matlab indexing!


				/*We are going to plug this vector into the earth model, at the right vertices corresponding to this particular 
				 * ice cap: */
				forcingglobal->SetValues(M,index,forcingfromcap,INS_VAL);
				xDelete<int>(index);
			}
		}

		/*Assemble vector:*/
		forcingglobal->Assemble();
		
		/*Plug into elements:*/
		InputUpdateFromVectorx(femmodel,forcingglobal,forcingenum,VertexSIdEnum);
	} 
	/*}}}*/

	/*Free ressources:{{{*/
	if(forcings){
		for(int i=0;i<nummodels-1;i++){
			IssmDouble* temp=forcings[i]; xDelete<IssmDouble>(temp);
		}
		xDelete<IssmDouble*>(forcings);
	}
	if(forcing)xDelete<IssmDouble>(forcing);
	if(forcingglobal)delete forcingglobal;
	if(transitions){
		for(int i=0;i<earthid;i++){
			IssmDouble* temp=transitions[i];
			xDelete<IssmDouble>(temp);
		}
		xDelete<IssmDouble*>(transitions);
		xDelete<int>(transitions_m);
		xDelete<int>(transitions_n);
	}
	if(nvs)xDelete<int>(nvs);
	/*}}}*/

} /*}}}*/
void TransferSealevel(FemModel* femmodel,int forcingenum){ /*{{{*/

	/*forcing being transferred from earth to ice caps: */
	IssmDouble*  forcing=NULL; 
	IssmDouble*  forcingglobal=NULL; 
	
	/*transition vectors:*/
	IssmDouble** transitions=NULL;
	int          ntransitions; 
	int*         transitions_m=NULL;
	int*         transitions_n=NULL;
	int          nv;
	
	/*communicators:*/
	ISSM_MPI_Comm fromcomm;
	ISSM_MPI_Comm* tocomms=NULL;
	ISSM_MPI_Status status;
	int         my_rank;
	int         modelid,earthid;
	int         nummodels;
	int         numcoms;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&modelid,ModelIdEnum);
	femmodel->parameters->FindParam(&earthid,EarthIdEnum);
	femmodel->parameters->FindParam(&nummodels,NumModelsEnum);
	my_rank=IssmComm::GetRank();
	
	/*retrieve the inter communicators that will be used to send data from earth to ice caps:*/
	if(modelid==earthid){
		GenericParam<ISSM_MPI_Comm*>* parcoms = dynamic_cast<GenericParam<ISSM_MPI_Comm*>*>(femmodel->parameters->FindParamObject(IcecapToEarthCommEnum));
		if(!parcoms)_error_("TransferSealevel error message: could not find IcecapToEarthComm communicator");
		tocomms=parcoms->GetParameterValue();
		//femmodel->parameters->FindParam((int**)(&tocomms),&numcoms,IcecapToEarthCommEnum);
	}
	else{
		GenericParam<ISSM_MPI_Comm>* parcom = dynamic_cast<GenericParam<ISSM_MPI_Comm>*>(femmodel->parameters->FindParamObject(IcecapToEarthCommEnum));
		if(!parcom)_error_("TransferSealevel error message: could not find IcecapToEarthComm communicator");
		fromcomm=parcom->GetParameterValue();
		//femmodel->parameters->FindParam((int*)(&fromcomm), IcecapToEarthCommEnum);
	}


	/*Retrieve sea-level on earth model: */
	if(modelid==earthid){
		nv=femmodel->vertices->NumberOfVertices();
		GetVectorFromInputsx(&forcingglobal,femmodel,forcingenum,VertexSIdEnum);
	}

	/*Send the forcing to the ice caps:{{{*/
	if(my_rank==0){
		
		if(modelid==earthid){
			
			/*Retrieve transition vectors, used to figure out global forcing contribution to each ice cap's own elements: */
			femmodel->parameters->FindParam(&transitions,&ntransitions,&transitions_m,&transitions_n,SealevelriseTransitionsEnum);
			
			if(ntransitions!=earthid)_error_("TransferSeaLevel error message: number of transition vectors is not equal to the number of icecaps!");

			for(int i=0;i<earthid;i++){
				nv=transitions_m[i];
				forcing=xNew<IssmDouble>(nv);
				IssmDouble* transition=transitions[i];
				for(int j=0;j<nv;j++){
					forcing[j]=forcingglobal[reCast<int>(transition[j])-1];
				}
				ISSM_MPI_Send(&nv, 1, ISSM_MPI_INT, 0, i, tocomms[i]);
				ISSM_MPI_Send(forcing, nv, ISSM_MPI_DOUBLE, 0, i, tocomms[i]);
			}
		}
		else{
			ISSM_MPI_Recv(&nv, 1, ISSM_MPI_INT, 0, modelid, fromcomm, &status);
			forcing=xNew<IssmDouble>(nv);
			ISSM_MPI_Recv(forcing, nv, ISSM_MPI_DOUBLE, 0, modelid, fromcomm, &status);
		}
	}
	/*}}}*/

	/*On each ice cap, spread the forcing across cpus, and update the elements dataset accordingly: {{{*/
	if(modelid!=earthid){

		ISSM_MPI_Bcast(&nv,1,ISSM_MPI_INT,0,IssmComm::GetComm());
		if(my_rank!=0)forcing=xNew<IssmDouble>(nv);
		ISSM_MPI_Bcast(forcing,nv,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

		/*Plug into elements:*/
		InputUpdateFromVectorx(femmodel,forcing,forcingenum,VertexSIdEnum);
	} 
	/*}}}*/

	/*Free ressources:{{{*/
	if(forcingglobal)xDelete<IssmDouble>(forcingglobal);
	if(forcing)xDelete<IssmDouble>(forcing);
	if(transitions){
		for(int i=0;i<ntransitions;i++){
			IssmDouble* temp=transitions[i];
			xDelete<IssmDouble>(temp);
		}
		xDelete<IssmDouble*>(transitions);
		xDelete<int>(transitions_m);
		xDelete<int>(transitions_n);
	}
	/*}}}*/

} /*}}}*/
