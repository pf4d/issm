/*!\file: CreateParametersOutputDefinitions.cpp
 * \brief driver for creating output definitions dataset, and including it into the parameters dataset
 */ 

#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void CreateOutputDefinitions(Elements* elements, Parameters* parameters,IoModel* iomodel){

	int i,j;
	
	DataSet*     output_definitions      = NULL;
	int*         output_definition_enums = NULL;
	int          num_output_definitions;

	/*Create output_definitions dataset: */
	output_definitions=new DataSet();

	char** out_strings = NULL;
	iomodel->FetchData(&out_strings,&num_output_definitions,"md.outputdefinition.list");
	if(num_output_definitions>0){
		output_definition_enums=xNew<int>(num_output_definitions);
		for(int i=0;i<num_output_definitions;i++){
			output_definition_enums[i]=StringToEnumx(out_strings[i]);
		}
	}
	// free data
	for(int i=0;i<num_output_definitions;i++) xDelete<char>(out_strings[i]);
	xDelete<char*>(out_strings);

	if(num_output_definitions){
		for (i=0;i<num_output_definitions;i++){
			if (output_definition_enums[i]==MassfluxatgateEnum){
				/*Deal with mass flux gates:{{{ */

				/*massfluxatgate variables: */
				int          temp,numgates;
				char       **gatenames           = NULL;
				char		  **gatedefinitionstrings = NULL;
				IssmDouble **gatesegments        = NULL;
				int         *gatesegments_M      = NULL;

				/*Fetch segments and names: */
				iomodel->FetchMultipleData(&gatenames,&numgates,                     "md.massfluxatgate.name");
				iomodel->FetchMultipleData(&gatedefinitionstrings,&temp,             "md.massfluxatgate.definitionstring"); _assert_(temp==numgates);
				iomodel->FetchMultipleData(&gatesegments,&gatesegments_M,NULL,&temp, "md.massfluxatgate.segments");       _assert_(temp==numgates);

				for(j=0;j<numgates;j++){
					output_definitions->AddObject(new Massfluxatgate<IssmDouble>(gatenames[j],StringToEnumx(gatedefinitionstrings[j]),gatesegments_M[j],gatesegments[j]));
				}
				/*Free ressources:*/
				for(j=0;j<numgates;j++){
					char*       string  = gatenames[j];             xDelete<char>(string);
					char*       string2 = gatedefinitionstrings[j]; xDelete<char>(string2);
					IssmDouble* gate    = gatesegments[j];          xDelete<IssmDouble>(gate);
				}
				xDelete<char*>(gatenames);
				xDelete<IssmDouble*>(gatesegments);
				xDelete<int>(gatesegments_M);
				xDelete<char*>(gatedefinitionstrings);
				/*}}}*/
			}
			else if (output_definition_enums[i]==MisfitEnum){
				/*Deal with misfits: {{{*/
				
				/*misfit variables: */
				int          nummisfits;
				char**       misfit_name_s						= NULL;    
				char**		 misfit_definitionstring_s		= NULL;    
				char**       misfit_model_string_s			= NULL;
				IssmDouble** misfit_observation_s			= NULL;
				char**		 misfit_observation_string_s	= NULL;
				int*         misfit_observation_M_s			= NULL;
				int*         misfit_observation_N_s			= NULL;
				int*         misfit_local_s					= NULL;
				char**       misfit_timeinterpolation_s	= NULL;
				IssmDouble** misfit_weights_s					= NULL;
				int*         misfit_weights_M_s				= NULL;
				int*         misfit_weights_N_s				= NULL;
				char**       misfit_weights_string_s		= NULL;

				/*Fetch name, model_string, observation, observation_string, etc ... (see src/m/classes/misfit.m): */
				iomodel->FetchMultipleData(&misfit_name_s,&nummisfits,                                                        "md.misfit.name");
				iomodel->FetchMultipleData(&misfit_definitionstring_s,&nummisfits,                                            "md.misfit.definitionstring");
				iomodel->FetchMultipleData(&misfit_model_string_s,&nummisfits,                                                "md.misfit.model_string");
				iomodel->FetchMultipleData(&misfit_observation_s,&misfit_observation_M_s,&misfit_observation_N_s,&nummisfits, "md.misfit.observation");
				iomodel->FetchMultipleData(&misfit_observation_string_s,&nummisfits,                                          "md.misfit.observation_string");
				iomodel->FetchMultipleData(&misfit_timeinterpolation_s,&nummisfits,                                           "md.misfit.timeinterpolation");
				iomodel->FetchMultipleData(&misfit_local_s,&nummisfits,                                                       "md.misfit.local");
				iomodel->FetchMultipleData(&misfit_weights_s,&misfit_weights_M_s,&misfit_weights_N_s,&nummisfits,             "md.misfit.weights");
				iomodel->FetchMultipleData(&misfit_weights_string_s,&nummisfits,                                              "md.misfit.weights_string");

				for(j=0;j<nummisfits;j++){

					/*First create a misfit object for that specific string (misfit_model_string_s[j]):*/
					output_definitions->AddObject(new Misfit(misfit_name_s[j],StringToEnumx(misfit_definitionstring_s[j]),StringToEnumx(misfit_model_string_s[j]),StringToEnumx(misfit_observation_string_s[j]),misfit_timeinterpolation_s[j],(bool)misfit_local_s[j],StringToEnumx(misfit_weights_string_s[j])));

					/*Now, for this particular misfit object, make sure we plug into the elements: the observation, and the weights.*/
					for(int k=0;k<elements->Size();k++){
						Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(k));
						element->InputCreate(misfit_observation_s[j], iomodel,misfit_observation_M_s[j],misfit_observation_N_s[j],1,StringToEnumx(misfit_observation_string_s[j]),7);
						element->InputCreate(misfit_weights_s[j], iomodel,misfit_weights_M_s[j],misfit_weights_N_s[j],1,StringToEnumx(misfit_weights_string_s[j]),7);
					}

				}

				/*Free ressources:*/
				for(j=0;j<nummisfits;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;

					string = misfit_definitionstring_s[j];		xDelete<char>(string);
					string = misfit_observation_string_s[j];	xDelete<char>(string);
					string = misfit_model_string_s[j];			xDelete<char>(string);
					string = misfit_weights_string_s[j];		xDelete<char>(string);
					string = misfit_name_s[j];    xDelete<char>(string);
					string = misfit_timeinterpolation_s[j];    xDelete<char>(string);
					matrix = misfit_observation_s[j]; xDelete<IssmDouble>(matrix);
					matrix = misfit_weights_s[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(misfit_name_s);
				xDelete<char*>(misfit_model_string_s);
				xDelete<char*>(misfit_definitionstring_s);
				xDelete<IssmDouble*>(misfit_observation_s);
				xDelete<char*>(misfit_observation_string_s);
				xDelete<int>(misfit_observation_M_s);
				xDelete<int>(misfit_observation_N_s);
				xDelete<int>(misfit_local_s);
				xDelete<char*>(misfit_timeinterpolation_s);
				xDelete<IssmDouble*>(misfit_weights_s);
				xDelete<int>(misfit_weights_M_s);
				xDelete<int>(misfit_weights_N_s);
				xDelete<char*>(misfit_weights_string_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==NodalvalueEnum){
				/*Deal with nodal values: {{{*/
				
				/*nodal value variables: */
				int          numnodalvalues;
				char**       nodalvalue_name_s             = NULL;    
				char**       nodalvalue_definitionstrings             = NULL;    
				char**       nodalvalue_modelstrings        = NULL;
				int*         nodalvalue_node_s = NULL;

				/*Fetch name, model_enum, etc ... (see src/m/classes/nodalvalue.m): */
				iomodel->FetchMultipleData(&nodalvalue_name_s,&numnodalvalues,            "md.nodalvalue.name");
				iomodel->FetchMultipleData(&nodalvalue_definitionstrings,&numnodalvalues, "md.nodalvalue.definitionenum");
				iomodel->FetchMultipleData(&nodalvalue_modelstrings,&numnodalvalues,      "md.nodalvalue.model_enum");
				iomodel->FetchMultipleData(&nodalvalue_node_s,&numnodalvalues,            "md.nodalvalue.node");

				for(j=0;j<numnodalvalues;j++){

					/*First create a nodalvalue object for that specific enum (nodalvalue_model_enum_s[j]):*/
					output_definitions->AddObject(new Nodalvalue(nodalvalue_name_s[j],StringToEnumx(nodalvalue_definitionstrings[j]),StringToEnumx(nodalvalue_modelstrings[j]),nodalvalue_node_s[j]-1)); //-1 because matlab to c indexing.
				}
					
				/*Free ressources:*/
				for(j=0;j<numnodalvalues;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;
					string = nodalvalue_name_s[j];    xDelete<char>(string);
				}
				xDelete<char*>(nodalvalue_name_s);
				xDelete<char*>(nodalvalue_modelstrings);
				xDelete<char*>(nodalvalue_definitionstrings);
				xDelete<int>(nodalvalue_node_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==MassconEnum){
				/*Deal with masscons: {{{*/

				/*masscon variables: */
				int          nummasscons;
				char**       masscon_name_s					= NULL;    
				int*         masscon_definitionenum_s		= NULL;    
				IssmDouble** masscon_levelset_s           = NULL;
				int*         masscon_levelset_M_s			= NULL;
				int*         masscon_levelset_N_s			= NULL;

				/*Fetch name and levelset, etc ... (see src/m/classes/masscon.m): */
				iomodel->FetchMultipleData(&masscon_name_s,&nummasscons,                                                "md.masscon.name");
				iomodel->FetchMultipleData(&masscon_definitionenum_s,&nummasscons,                                      "md.masscon.definitionenum");
				iomodel->FetchMultipleData(&masscon_levelset_s,&masscon_levelset_M_s,&masscon_levelset_N_s,&nummasscons,"md.masscon.levelset");
				for(j=0;j<nummasscons;j++){

					/*Create a masscon object: */
					output_definitions->AddObject(new Masscon(masscon_name_s[j],masscon_definitionenum_s[j],masscon_levelset_s[j],masscon_levelset_M_s[j]));

				}

				/*Free ressources:*/
				for(j=0;j<nummasscons;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;

					string = masscon_name_s[j];    xDelete<char>(string);
					matrix = masscon_levelset_s[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(masscon_name_s);
				xDelete<IssmDouble*>(masscon_levelset_s);
				xDelete<int>(masscon_levelset_M_s);
				xDelete<int>(masscon_levelset_N_s);
				xDelete<int>(masscon_definitionenum_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==MassconaxpbyEnum){
				/*Deal with masscon combinations: {{{*/
				
				/*masscon variables: */
				char**       masscon_name_s             = NULL;    
				int*         masscon_definitionenum_s             = NULL;    
				char**       masscon_namex_s             = NULL;    
				char**       masscon_namey_s             = NULL;    
				IssmDouble*  masscon_alpha_s     = NULL;
				IssmDouble*  masscon_beta_s     = NULL;
				int          num;

				/*Fetch names and multiplicators, etc ... (see src/m/classes/masscon_axpby.m): */
				iomodel->FetchMultipleData(&masscon_name_s,&num,          "md.massconaxpby.name");
				iomodel->FetchMultipleData(&masscon_definitionenum_s,&num,"md.massconaxpby.definitionenum");
				iomodel->FetchMultipleData(&masscon_namex_s,&num,         "md.massconaxpby.namex");
				iomodel->FetchMultipleData(&masscon_namey_s,&num,         "md.massconaxpby.namey");
				iomodel->FetchMultipleData(&masscon_alpha_s,&num,         "md.massconaxpby.alpha");
				iomodel->FetchMultipleData(&masscon_beta_s,&num,          "md.massconaxpby.beta");
				for(j=0;j<num;j++){

					/*Create a masscon axpyb object: */
					output_definitions->AddObject(new Massconaxpby(masscon_name_s[j],masscon_definitionenum_s[j],masscon_namex_s[j],masscon_namey_s[j],masscon_alpha_s[j],masscon_beta_s[j]));

				}

				/*Free ressources:*/
				for(j=0;j<num;j++){
					char* string=NULL;
					string = masscon_name_s[j];    xDelete<char>(string);
					string = masscon_namex_s[j];    xDelete<char>(string);
					string = masscon_namey_s[j];    xDelete<char>(string);
				}
				xDelete<char*>(masscon_name_s);
				xDelete<char*>(masscon_namex_s);
				xDelete<char*>(masscon_namey_s);
				xDelete<int>(masscon_definitionenum_s);
				xDelete<IssmDouble>(masscon_alpha_s);
				xDelete<IssmDouble>(masscon_beta_s);
				/*}}}*/
			}
			else _error_("output definition enum " << output_definition_enums[i] << " not supported yet!");
		}
	}
	parameters->AddObject(new DataSetParam(OutputdefinitionEnum,output_definitions));

	/*Free ressources:*/
	delete output_definitions;
	xDelete<int>(output_definition_enums);
}
