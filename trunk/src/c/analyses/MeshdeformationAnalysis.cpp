#include "./MeshdeformationAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void MeshdeformationAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void MeshdeformationAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void MeshdeformationAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
int  MeshdeformationAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void MeshdeformationAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void MeshdeformationAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/

/*Finite Element Analysis*/
void           MeshdeformationAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* MeshdeformationAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* MeshdeformationAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* MeshdeformationAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* MeshdeformationAnalysis::CreatePVector(Element* element){/*{{{*/
_error_("not implemented yet");
}/*}}}*/
void           MeshdeformationAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           MeshdeformationAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           MeshdeformationAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           MeshdeformationAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
