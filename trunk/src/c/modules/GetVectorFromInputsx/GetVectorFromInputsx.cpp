/*!\file GetVectorFromInputsx
 * \brief retrieve vector from inputs in elements
 */

#include "./GetVectorFromInputsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void GetVectorFromInputsx(Vector<IssmDouble>** pvector,FemModel* femmodel,int name,int type){

	int i;
	Vector<IssmDouble>* vector=NULL;

	switch(type){
		case ElementSIdEnum:
			vector=new Vector<IssmDouble>(femmodel->elements->NumberOfElements());
			break;
		case VertexPIdEnum: case VertexSIdEnum:
			vector=new Vector<IssmDouble>(femmodel->vertices->NumberOfVertices());
			break;
		case NodesEnum:case NodeSIdEnum:
			vector=new Vector<IssmDouble>(femmodel->nodes->NumberOfNodes());
			break;
		default:
			_error_("vector type: " << EnumToStringx(type) << " not supported yet!");
	}
	/*Look up in elements*/
	for(i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->GetVectorFromInputs(vector,name,type);
	}

	vector->Assemble();

	/*Assign output pointers:*/
	*pvector=vector;
}

void GetVectorFromInputsx(IssmDouble** pvector,FemModel* femmodel,int name, int type){

	/*output: */
	IssmDouble* vector=NULL;

	/*intermediary: */
	Vector<IssmDouble>* vec_vector=NULL;

	GetVectorFromInputsx(&vec_vector,femmodel,name,type);
	vector=vec_vector->ToMPISerial();

	/*Free ressources:*/
	delete vec_vector;

	/*Assign output pointers:*/
	*pvector=vector;
}
