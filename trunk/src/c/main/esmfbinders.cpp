/*!\file:  esmfbinder.cpp
 * \brief: ESMF binders for ISSM. Binders developed initially for the GEOS-5 framework.
 */ 


#include "./issm.h"

/*GCM specific declarations:*/
const int GCMForcingNumTerms = 1;
const int GCMForcingTerms[GCMForcingNumTerms]= { SMBgcmEnum}; 

const int ISSMOutputNumTerms = 1;
const int ISSMOutputTerms[ISSMOutputNumTerms]= { SurfaceEnum };

extern "C" {

	FemModel *femmodel;

	void InitializeISSM(int argc, char** argv, int** pelementsonlocalrank, int* pnumberofelements, ISSM_MPI_Comm comm_init){ /*{{{*/

		int numberofelements;
		int* elementsonlocalrank=NULL;

		/*Initialize femmodel from arguments provided command line: */
		femmodel = new FemModel(argc,argv,comm_init);

		/*Figure out the partition for elements, and return: */
		numberofelements=femmodel->elements->NumberOfElements();

		elementsonlocalrank=xNewZeroInit<int>(numberofelements); 
		for (int i=0;i<femmodel->elements->Size();i++){
			Element* element=dynamic_cast<Element*>(femmodel->elements->GetObjectByOffset(i));
			elementsonlocalrank[element->Sid()]=1;
		}

		/*Some specific code here for the binding: */
		femmodel->parameters->SetParam(SMBgcmEnum,SmbEnum); //bypass SMB model, will be provided by GCM!
	
		/*Restart file: */
		femmodel->Restart();

		/*Assign output pointers: */
		*pnumberofelements=numberofelements;
		*pelementsonlocalrank=elementsonlocalrank;

	} /*}}}*/
	void RunISSM(IssmDouble dt, IssmDouble* gcmforcings, IssmDouble* issmoutputs){ /*{{{*/

		int numberofelements;
		IssmDouble yts;
		IssmDouble rho_ice;
		IssmDouble area;
		IssmDouble start_time,final_time;

		/*Figure out number of elements: */
		numberofelements=femmodel->elements->Size();

		/*Fetch some necessary constants: */
		femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);

		/*Setup gcm forcings as element-wise input: {{{ */
		for (int f=0;f<GCMForcingNumTerms;f++){

			int forcing_type=GCMForcingTerms[f];

			for (int i=0;i<femmodel->elements->Size();i++){
				Element* element=dynamic_cast<Element*>(femmodel->elements->GetObjectByOffset(i));

				switch(forcing_type){
					case SMBgcmEnum:
						/*{{{*/
						{

						/*Recover rho_ice: */
						rho_ice=element->matpar->GetMaterialParameter(MaterialsRhoIceEnum);

						/*Recover area of element: */
						area=element->SurfaceArea();

						/*Recover smb forcing from the gcm forcings: */
						IssmDouble smbforcing=*(gcmforcings+f*numberofelements+i); 

						/*Convert to SI. The smbforcing from GEOS-5 in kg/s, and we transform it into m/s: */
						smbforcing=smbforcing/(rho_ice*area);

						/*Add into the element as new forcing :*/
						element->inputs->AddInput(new DoubleInput(SmbMassBalanceEnum,smbforcing));

						}
						/*}}}*/
						break; 
					default: 
						{ _error_("Unknown forcing type " << forcing_type << "\n"); }
						break;
				}
			}
		}

		/*}}}*/

		/*Retrieve ISSM outputs and pass them back to the Gcm : {{{*/
		for (int f=0;f<ISSMOutputNumTerms;f++){

			int output_type=ISSMOutputTerms[f];

			for (int i=0;i<femmodel->elements->Size();i++){
				Element* element=dynamic_cast<Element*>(femmodel->elements->GetObjectByOffset(i));

				switch(output_type){
					case SurfaceEnum:
						/*{{{*/
						{

						IssmDouble surface;
						
						/*Recover surface from the ISSM element: */
						Input* surface_input = element->GetInput(SurfaceEnum); _assert_(surface_input);
						surface_input->GetInputAverage(&surface);
			
						*(issmoutputs+f*numberofelements+i) = surface;

						}
						/*}}}*/
						break; 
					default: 
						{ _error_("Unknown output type " << output_type << "\n"); }
						break;
				}
			}
		}

		/*}}}*/

		/*Before running, setup the time interval: */
		femmodel->parameters->FindParam(&start_time,TimeEnum);
		final_time=start_time+dt;
		femmodel->parameters->SetParam(final_time,TimesteppingFinalTimeEnum); //we are bypassing ISSM's initial final time!

		/*Now, run: */
		femmodel->Solve();

		/*For the next time around, save the final time as start time */
		femmodel->parameters->SetParam(final_time,TimeEnum);
	} /*}}}*/
	void FinalizeISSM(){ /*{{{*/

		/*Output results: */
		OutputResultsx(femmodel);
			
		/*Check point: */
		femmodel->CheckPoint();

		/*Wrap up: */
		delete femmodel; femmodel=NULL;
	} /*}}}*/

}
