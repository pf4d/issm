#include "./SealevelriseAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void SealevelriseAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*No constraints*/
}/*}}}*/
void SealevelriseAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*No loads*/
}/*}}}*/
void SealevelriseAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/
	::CreateNodes(nodes,iomodel,SealevelriseAnalysisEnum,P1Enum);
}/*}}}*/
int  SealevelriseAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void SealevelriseAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	/*Create inputs: */
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.mask.land_levelset",MaskLandLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.slr.deltathickness",SealevelriseDeltathicknessEnum);
	iomodel->FetchDataToInput(elements,"md.slr.sealevel",SealevelEnum,0);

}/*}}}*/
void SealevelriseAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int         nl;
	IssmDouble* love_h=NULL;
	IssmDouble* love_k=NULL;
	IssmDouble* love_l=NULL;
	
	bool elastic=false;
	IssmDouble* G_elastic = NULL;
	IssmDouble* G_elastic_local = NULL;
	IssmDouble* U_elastic = NULL;
	IssmDouble* U_elastic_local = NULL;
	IssmDouble* H_elastic = NULL;
	IssmDouble* H_elastic_local = NULL;
	int         M,m,lower_row,upper_row;
	IssmDouble  degacc=.01;

	int     numoutputs;
	char**  requestedoutputs = NULL;

	/*transition vectors: */
	IssmDouble **transitions    = NULL;
	int         *transitions_M    = NULL;
	int         *transitions_N    = NULL;
	int          ntransitions;

	/*some constant parameters: */
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.reltol",SealevelriseReltolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.abstol",SealevelriseAbstolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.maxiter",SealevelriseMaxiterEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.rigid",SealevelriseRigidEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.elastic",SealevelriseElasticEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.rotation",SealevelriseRotationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.tide_love_h",SealevelriseTidalLoveHEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.tide_love_k",SealevelriseTidalLoveKEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.fluid_love",SealevelriseFluidLoveEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.equatorial_moi",SealevelriseEquatorialMoiEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.polar_moi",SealevelrisePolarMoiEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.angular_velocity",SealevelriseAngularVelocityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.slr.ocean_area_scaling",SealevelriseOceanAreaScalingEnum));

	iomodel->FetchData(&elastic,"md.slr.elastic");
	if(elastic){

		/*love numbers: */
		iomodel->FetchData(&love_h,&nl,NULL,"md.slr.love_h");
		iomodel->FetchData(&love_k,&nl,NULL,"md.slr.love_k");
		iomodel->FetchData(&love_l,&nl,NULL,"md.slr.love_l");

		/*compute elastic green function for a range of angles*/
		iomodel->FetchData(&degacc,"md.slr.degacc");
		M=reCast<int,IssmDouble>(180./degacc+1.);
		G_elastic=xNew<IssmDouble>(M);
		U_elastic=xNew<IssmDouble>(M);
		H_elastic=xNew<IssmDouble>(M);
		
		/*compute combined legendre + love number (elastic green function:*/
		m=DetermineLocalSize(M,IssmComm::GetComm());
		GetOwnershipBoundariesFromRange(&lower_row,&upper_row,m,IssmComm::GetComm());
		G_elastic_local=xNew<IssmDouble>(m);
		U_elastic_local=xNew<IssmDouble>(m);
		H_elastic_local=xNew<IssmDouble>(m);

		for(int i=lower_row;i<upper_row;i++){
			IssmDouble alpha,x;
			alpha= reCast<IssmDouble>(i)*degacc * PI / 180.0;

			G_elastic_local[i-lower_row]= (love_k[nl-1]-love_h[nl-1])/2.0/sin(alpha/2.0);
			U_elastic_local[i-lower_row]= (love_h[nl-1])/2.0/sin(alpha/2.0);
			H_elastic_local[i-lower_row]= 0; 
			IssmDouble Pn,Pn1,Pn2;
			IssmDouble Pn_p,Pn_p1,Pn_p2;
			for (int n=0;n<nl;n++) {
				IssmDouble deltalove_G;
				IssmDouble deltalove_U;

				deltalove_G = (love_k[n]-love_k[nl-1]-love_h[n]+love_h[nl-1]);
				deltalove_U = (love_h[n]-love_h[nl-1]);
		
				/*compute legendre polynomials: P_n(cos\theta) & d P_n(cos\theta)/ d\theta: */
				if(n==0){
					Pn=1; 
					Pn_p=0; 
				}
				else if(n==1){ 
					Pn = cos(alpha); 
					Pn_p = 1; 
				}
				else{
					Pn = ( (2*n-1)*cos(alpha)*Pn1 - (n-1)*Pn2 ) /n;
					Pn_p = ( (2*n-1)*(Pn1+cos(alpha)*Pn_p1) - (n-1)*Pn_p2 ) /n;
				}
				Pn2=Pn1; Pn1=Pn;
				Pn_p2=Pn_p1; Pn_p1=Pn_p;

				G_elastic_local[i-lower_row] += deltalove_G*Pn;		// gravitational potential 
				U_elastic_local[i-lower_row] += deltalove_U*Pn;		// vertical (up) displacement 
				H_elastic_local[i-lower_row] += sin(alpha)*love_l[n]*Pn_p;		// horizontal displacements 
			}
		}

		/*merge G_elastic_local into G_elastic; U_elastic_local into U_elastic; H_elastic_local to H_elastic:{{{*/
		int* recvcounts=xNew<int>(IssmComm::GetSize());
		int* displs=xNew<int>(IssmComm::GetSize());

		//recvcounts:
		ISSM_MPI_Allgather(&m,1,ISSM_MPI_INT,recvcounts,1,ISSM_MPI_INT,IssmComm::GetComm());

		/*displs: */
		ISSM_MPI_Allgather(&lower_row,1,ISSM_MPI_INT,displs,1,ISSM_MPI_INT,IssmComm::GetComm());

		/*All gather:*/
		ISSM_MPI_Allgatherv(G_elastic_local, m, ISSM_MPI_DOUBLE, G_elastic, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
		ISSM_MPI_Allgatherv(U_elastic_local, m, ISSM_MPI_DOUBLE, U_elastic, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
		ISSM_MPI_Allgatherv(H_elastic_local, m, ISSM_MPI_DOUBLE, H_elastic, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
		/*free ressources: */
		xDelete<int>(recvcounts);
		xDelete<int>(displs);

		/*}}}*/

		/*Avoid singularity at 0: */
		G_elastic[0]=G_elastic[1];
		parameters->AddObject(new DoubleVecParam(SealevelriseGElasticEnum,G_elastic,M));
		U_elastic[0]=U_elastic[1];
		parameters->AddObject(new DoubleVecParam(SealevelriseUElasticEnum,U_elastic,M));
		H_elastic[0]=H_elastic[1];
		parameters->AddObject(new DoubleVecParam(SealevelriseHElasticEnum,H_elastic,M));

		/*free ressources: */
		xDelete<IssmDouble>(love_h);
		xDelete<IssmDouble>(love_k);
		xDelete<IssmDouble>(love_l);
		xDelete<IssmDouble>(G_elastic);
		xDelete<IssmDouble>(G_elastic_local);
		xDelete<IssmDouble>(U_elastic);
		xDelete<IssmDouble>(U_elastic_local);
		xDelete<IssmDouble>(H_elastic);
		xDelete<IssmDouble>(H_elastic_local);
	}
	
	/*Transitions: */
	iomodel->FetchData(&transitions,&transitions_M,&transitions_N,&ntransitions,"md.slr.transitions");
	if(transitions){
		parameters->AddObject(new DoubleMatArrayParam(SealevelriseTransitionsEnum,transitions,ntransitions,transitions_M,transitions_N));

		for(int i=0;i<ntransitions;i++){
			IssmDouble* transition=transitions[i];
			xDelete<IssmDouble>(transition);
		}
		xDelete<IssmDouble*>(transitions);
		xDelete<int>(transitions_M);
		xDelete<int>(transitions_N);
	}

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.slr.requested_outputs");
	parameters->AddObject(new IntParam(SealevelriseNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(SealevelriseRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.slr.requested_outputs");


}/*}}}*/

/*Finite Element Analysis*/
void           SealevelriseAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* SealevelriseAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* SealevelriseAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* SealevelriseAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* SealevelriseAnalysis::CreatePVector(Element* element){/*{{{*/
_error_("not implemented yet");
}/*}}}*/
void           SealevelriseAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           SealevelriseAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           SealevelriseAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	
	IssmDouble *deltaS  = NULL;
	IssmDouble *S  = NULL;
	int*        sidlist = NULL;
	int         numvertices;
	
	numvertices= element->GetNumberOfVertices();
	sidlist=xNew<int>(numvertices);
	
	element->GetVerticesSidList(sidlist);

	deltaS = xNew<IssmDouble>(numvertices);
	for(int i=0;i<numvertices;i++){
		deltaS[i]=solution[sidlist[i]];
	}

	S = xNew<IssmDouble>(numvertices);
	element->GetInputListOnVertices(S,SealevelEnum,0);

	/*Add deltaS to S:*/
	for (int i=0;i<numvertices;i++)S[i]+=deltaS[i];

	/*Add S back into inputs: */
	element->AddInput(SealevelEnum,S,P1Enum);

	/*Free ressources:*/
	xDelete<int>(sidlist);
	xDelete<IssmDouble>(deltaS);
	xDelete<IssmDouble>(S);

}/*}}}*/
void           SealevelriseAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
