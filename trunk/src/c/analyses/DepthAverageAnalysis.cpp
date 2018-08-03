#include "./DepthAverageAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void DepthAverageAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
}/*}}}*/
void DepthAverageAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
}/*}}}*/
void DepthAverageAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	::CreateNodes(nodes,iomodel,DepthAverageAnalysisEnum,P1Enum);

}/*}}}*/
int  DepthAverageAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void DepthAverageAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	if(iomodel->domaintype==Domain2DverticalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
	}
}/*}}}*/
void DepthAverageAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           DepthAverageAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* DepthAverageAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* DepthAverageAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* DepthAverageAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries */
	int         dim;
	IssmDouble  Jdet,D,dt=1.e+9;
	IssmDouble *xyz_list = NULL;

	/*Get dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    B      = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		GetB(B,element,dim,xyz_list,gauss);
		D=gauss->weight*Jdet*dt;

		/*vertical diffusion*/
		TripleMultiply(B,1,numnodes,1,
					&D,1,1,0,
					B,1,numnodes,0,
					&Ke->values[0],1);

		/*Next value*/
		D=gauss->weight*Jdet;
		element->NodalFunctions(B,gauss);
		TripleMultiply(B,numnodes,1,0,
					&D,1,1,0,
					B,1,numnodes,0,
					&Ke->values[0],1);
	} 

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(B);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* DepthAverageAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries*/
	int         input_enum;
	IssmDouble  Jdet,scalar,value;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe     = element->NewElementVector();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&input_enum,InputToDepthaverageInEnum);
	Input* input = element->GetInput(input_enum); _assert_(input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(3);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		/* Build transient now */
		input->GetInputValue(&value, gauss);
		scalar=value*Jdet*gauss->weight;
		for(int i=0;i<numnodes;i++) pe->values[i]+=scalar*basis[i];

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return pe;
}/*}}}*/
void           DepthAverageAnalysis::GetB(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*	Compute B  matrix. B=[dh1/dz dh2/dz dh3/dz dh4/dz dh5/dz dh6/dz];
		where hi is the interpolation function for node i.*/

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* dbasis=xNew<IssmDouble>(dim*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B: */
	for(int i=0;i<numnodes;i++){
		B[i] = dbasis[(dim-1)*numnodes+i];
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);
}
/*}}}*/
void           DepthAverageAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           DepthAverageAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           DepthAverageAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int inputenum;
	element->FindParam(&inputenum,InputToDepthaverageOutEnum);
	element->InputUpdateFromSolutionOneDof(solution,inputenum);
}/*}}}*/
void           DepthAverageAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
