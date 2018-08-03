/*!\file ConfigureObjectsx
 * \brief: configure objects in elements and loads to link in with nodes
 */

#include "./ConfigureObjectsx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"

int	ConfigureObjectsx( Elements* elements, Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters){

	/*Intermediary*/
	int       i;
	int       noerr = 1;
	int       configuration_type;
	Element  *element            = NULL;
	Load     *load               = NULL;
	Material *material           = NULL;

	/*Get analysis type: */
	parameters->FindParam(&configuration_type,ConfigurationTypeEnum);

	if(VerboseMProcessor()) _printf0_("      Configuring elements...\n");
	for(i=0;i<elements->Size();i++){
		element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->Configure(elements,loads,nodes,vertices,materials,parameters);
	}
	if(VerboseMProcessor()) _printf0_("      Configuring loads...\n");
	for(i=0;i<loads->Size();i++){
		load=(Load*)loads->GetObjectByOffset(i);
		if (load->InAnalysis(configuration_type)){
			load->Configure(elements,loads,nodes,vertices,materials,parameters);
		}
	}
	if(VerboseMProcessor()) _printf0_("      Configuring materials...\n");
	for(i=0;i<materials->Size();i++){
		material=(Material*)materials->GetObjectByOffset(i);
		material->Configure(elements);
	}
	return noerr;
}
