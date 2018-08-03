/*!\file OutputDefinitionsResponsex
 * \brief retrieve vector from inputs in elements
 */

#include "./OutputDefinitionsResponsex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"

IssmDouble OutputDefinitionsResponsex(FemModel* femmodel,const char* output_string){

	/*Ok, go find the output definitions dataset in the parameters, where our responses are hiding: */
	DataSet* output_definitions=((DataSetParam*)femmodel->parameters->FindParamObject(OutputdefinitionEnum))->value;

	/*Now, go through the output definitions, and retrieve the object which corresponds to our requested response, output_string: */
	for(int i=0;i<output_definitions->Size();i++){
		Definition* definition=dynamic_cast<Definition*>(output_definitions->GetObjectByOffset(i));

		char* name = definition->Name();
		if(strcmp(name,output_string)==0){

			/*This is the object that we have been chasing for. compute the response and return: */
			IssmDouble return_value=definition->Response(femmodel);
		
			/*cleanup: */
			xDelete<char>(name);

			/*return:*/
			return return_value;
		}
		xDelete<char>(name);
	}
	
	/*If we are here, did not find the definition for this response, not good!: */
	_error_("Could not find the response for output definition " << output_string << " because could not find the definition itself!");

}

IssmDouble OutputDefinitionsResponsex(FemModel* femmodel,int output_enum){

	/*Ok, go find the output definitions dataset in the parameters, where our responses are hiding: */
	DataSet* output_definitions=((DataSetParam*)femmodel->parameters->FindParamObject(OutputdefinitionEnum))->value;

	/*Now, go through the output definitions, and retrieve the object which corresponds to our requested response, output_enum: */
	for(int i=0;i<output_definitions->Size();i++){
		
		//Definition* definition=xDynamicCast<Definition*>(output_definitions->GetObjectByOffset(i));
		Definition* definition=dynamic_cast<Definition*>(output_definitions->GetObjectByOffset(i));

		int en = definition->DefinitionEnum();
		if(en==output_enum){

			/*This is the object that we have been chasing for. compute the response and return: */
			IssmDouble return_value=definition->Response(femmodel);
		
			/*return:*/
			return return_value;
		}
	}
	
	/*If we are here, did not find the definition for this response, not good!: */
	_error_("Could not find the response for output definition " << output_enum << " because could not find the definition itself!");

}
