/*!\file InputUpdateFromMatrixDakotax
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromMatrixDakotax.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../InputUpdateFromVectorDakotax/InputUpdateFromVectorDakotax.h"

void InputUpdateFromMatrixDakotax(FemModel* femmodel,double* matrix,int nrows,int ncols, int name, int type){

	int i;
	int numberofvertices;

	numberofvertices=femmodel->vertices->NumberOfVertices();

	if((ncols==1) && (nrows==numberofvertices)) InputUpdateFromVectorDakotax(femmodel,matrix,name,type);
	else{

		/*Update elements, nodes, loads and materials from inputs: */
		for(i=0;i<femmodel->elements->Size();i++){
			Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
			element->InputUpdateFromMatrixDakota(matrix,nrows,ncols,name,type);
		}
		for(i=0;i<femmodel->loads->Size();i++){
			Load* load=(Load*)femmodel->loads->GetObjectByOffset(i);
			load->InputUpdateFromMatrixDakota(matrix,nrows,ncols,name,type);
		}
		for(i=0;i<femmodel->materials->Size();i++){
			Material* material=(Material*)femmodel->materials->GetObjectByOffset(i);
			material->InputUpdateFromMatrixDakota(matrix,nrows,ncols,name,type);
		}
	}
}
