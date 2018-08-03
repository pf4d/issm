#include "./MeltingAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void MeltingAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*No Constraints*/
}/*}}}*/
void MeltingAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*if 2d: Error*/
	if(iomodel->domaintype==Domain2DhorizontalEnum) _error_("2d meshes not supported yet");

	//create penalties for nodes: no node can have a temperature over the melting point
	iomodel->FetchData(1,"md.mesh.vertexonbase");
	CreateSingleNodeToElementConnectivity(iomodel);

	for(int i=0;i<iomodel->numberofvertices;i++){
		if(iomodel->my_vertices[i]){
			if (reCast<int>(iomodel->Data("md.mesh.vertexonbase")[i])){
				loads->AddObject(new Pengrid(iomodel->loadcounter+i+1,i,iomodel,MeltingAnalysisEnum));
			}
		}
	}
	iomodel->DeleteData(1,"md.mesh.vertexonbase");

}/*}}}*/
void MeltingAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,MeltingAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  MeltingAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void MeltingAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int frictionlaw;

	/*Now, is the model 3d? otherwise, do nothing: */
	if(iomodel->domaintype==Domain2DhorizontalEnum)return;

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
	iomodel->FetchDataToInput(elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(elements,"md.slr.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(elements,"md.initialization.pressure",PressureEnum);
}/*}}}*/
void MeltingAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           MeltingAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* MeltingAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* MeltingAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* MeltingAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Get basal element*/
	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	IssmDouble  D,Jdet;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke    = basalelement->NewElementMatrix(NoneApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	IssmDouble latentheat   = element->GetMaterialParameter(MaterialsLatentheatEnum);
	IssmDouble heatcapacity = element->GetMaterialParameter(MaterialsHeatcapacityEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);
		D=latentheat/heatcapacity*gauss->weight*Jdet;

		TripleMultiply(basis,1,numnodes,1,
					&D,1,1,0,
					basis,1,numnodes,0,
					&Ke->values[0],1);
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	basalelement->DeleteMaterials(); delete basalelement;
	return Ke;
}/*}}}*/
ElementVector* MeltingAnalysis::CreatePVector(Element* element){/*{{{*/
	return NULL;
}/*}}}*/
void           MeltingAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           MeltingAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           MeltingAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	element->InputUpdateFromSolutionOneDof(solution,BasalforcingsGroundediceMeltingRateEnum);
}/*}}}*/
void           MeltingAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
