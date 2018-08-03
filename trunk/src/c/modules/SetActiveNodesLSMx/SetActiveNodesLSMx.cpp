/*!\file GetMaskOfIceVerticesLSMx 
 * \brief: Return a mask for all the vertices determining whether the node should be active or not. 
 */

#include "./SetActiveNodesLSMx.h"

#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../modules.h"

void SetActiveNodesLSMx(FemModel* femmodel){/*{{{*/
	/* activate/deactivate nodes for levelset method according to IceMaskNodeActivation */

	/* intermediaries */
	bool solvein2d=false;
	int i,in,domaintype,analysis_type;
	Elements* elements = femmodel->elements;

	/* find parameters */
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);

	for(i=0;i<elements->Size();i++){
		Element    *element  = xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		int         numnodes = element->GetNumberOfNodes();
		IssmDouble *mask     = xNew<IssmDouble>(numnodes);
		// include switch for elements with multiple different sets of nodes
		switch(element->GetElementType()){
			case MINIEnum:case MINIcondensedEnum:
			case TaylorHoodEnum:case XTaylorHoodEnum:case LATaylorHoodEnum:
			case CrouzeixRaviartEnum:case LACrouzeixRaviartEnum:case OneLayerP4zEnum:{
				Input* input=element->GetInput(IceMaskNodeActivationEnum);
				if(!input) _error_("Input " << EnumToStringx(IceMaskNodeActivationEnum) << " not found in element");

				/* Start looping on the number of vertices: */
				Gauss* gauss=element->NewGauss();
				for(int iv=0;iv<element->NumberofNodesVelocity();iv++){
					gauss->GaussNode(element->VelocityInterpolation(),iv);
					input->GetInputValue(&mask[iv],gauss);
				}
				for(int iv=0;iv<element->NumberofNodesPressure();iv++){
					gauss->GaussNode(element->PressureInterpolation(),iv);
					input->GetInputValue(&mask[element->NumberofNodesVelocity()+iv],gauss);
				}
				delete gauss;
				break;
			}
			default:
				element->GetInputListOnNodes(&mask[0],IceMaskNodeActivationEnum);
				break;
		}

		for(in=0;in<numnodes;in++){
			Node* node=element->GetNode(in);
			if(mask[in]==1.){
				node->Activate();
			}
			else {
				node->Deactivate();
			}
		}

		xDelete<IssmDouble>(mask);
	}
}/*}}}*/
void GetMaskOfIceVerticesLSMx(FemModel* femmodel){/*{{{*/

	/* Intermediaries */
	int i;

	/*Initialize vector with number of vertices*/
	int numvertices=femmodel->vertices->NumberOfVertices();
	Vector<IssmDouble>* vec_mask_ice=new Vector<IssmDouble>(numvertices); //vertices that have ice at next time step
	/*Fill vector with values: */
	for(i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		SetMaskOfIceElement(vec_mask_ice, element);
	}

	/*Assemble vector and serialize */
	vec_mask_ice->Assemble();
	IssmDouble* mask_ice=vec_mask_ice->ToMPISerial();
	InputUpdateFromVectorx(femmodel,mask_ice,IceMaskNodeActivationEnum,VertexSIdEnum);

	/*Clean up and return*/
	delete vec_mask_ice;
	xDelete<IssmDouble>(mask_ice);

}/*}}}*/
void SetMaskOfIceElement(Vector<IssmDouble>* vec_mask_ice, Element* element){/*{{{*/

	/* Intermediaries */
	int numvertices = element->GetNumberOfVertices();
	
	if(element->IsIceInElement()){
		for(int i = 0;i<numvertices;i++){
			vec_mask_ice->SetValue(element->vertices[i]->Sid(),1.,INS_VAL);
		}
	}
}/*}}}*/
