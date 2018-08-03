/*!\file GetVectorFromControlInputsx
 * \brief retrieve vector from inputs in elements
 */

#include "./GetVectorFromControlInputsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void GetVectorFromControlInputsx(Vector<IssmDouble>** pvector, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters,const char* data,bool onsid){

	int  num_controls;
	int *control_type = NULL;
	Vector<IssmDouble>*  vector=NULL;

	/*Retrieve some parameters*/
	parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);

	/*Allocate and populate gradient*/
	vector=new Vector<IssmDouble>(num_controls*vertices->NumberOfVertices());

	for(int i=0;i<num_controls;i++){
		for(int j=0;j<elements->Size();j++){
			Element* element=(Element*)elements->GetObjectByOffset(j);
			element->GetVectorFromControlInputs(vector,control_type[i],i,data,onsid);
		}
	}

	vector->Assemble();

	/*Assign output pointers:*/
	xDelete<int>(control_type);
	*pvector=vector;
}

void GetVectorFromControlInputsx( IssmDouble** pvector, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters, const char* data,bool onsid){

	/*output: */
	IssmDouble* vector=NULL;

	/*intermediary: */
	Vector<IssmDouble>* vec_vector=NULL;

	GetVectorFromControlInputsx( &vec_vector, elements,nodes, vertices, loads, materials, parameters,data,onsid);
	vector=vec_vector->ToMPISerial();

	/*Free ressources:*/
	delete vec_vector;

	/*Assign output pointers:*/
	*pvector=vector;
}
