/* \file ModelProcessorx.h
 * \brief  Header file for model processor
 */

#ifndef _MODEL_PROCESSORX_H_
#define _MODEL_PROCESSORX_H_

#include "../../classes/classes.h"
#include "../../analyses/analyses.h"

void ModelProcessorx(Elements** pelements, Nodes** pnodes, Vertices** pvertices, Materials** pmaterials, Constraints** pconstraints, Loads** ploads, Parameters** pparameters,IoModel* iomodel,FILE* toolkitfile, char* rootpath,const int solution_type,const int nummodels,const int* analysis_type_listh);

/*Creation of fem datasets: general drivers*/
void CreateElementsVerticesAndMaterials(Elements* elements,Vertices* vertices,Materials* materials, IoModel* iomodel,const int nummodels);
void CreateParameters(Parameters*parameters,IoModel* iomodel,char* rootpath,FILE* toolkitfile,const int solution_type);
void CreateParametersAutodiff(Parameters* parameters,IoModel* iomodel);
void CreateParametersControl(Parameters* parameters,IoModel* iomodel,int solution_type);
void CreateParametersDakota(Parameters* parameters,IoModel* iomodel,char* rootpath);
void CreateOutputDefinitions(Elements* elements, Parameters* parameters,IoModel* iomodel);
void UpdateElementsAndMaterialsControl(Elements* elements,Materials* materials, IoModel* iomodel);
void UpdateElementsAndMaterialsDakota(Elements* elements,Materials* materials, IoModel* iomodel);
void UpdateElementsTransient(Elements* elements,Parameters* parameters,IoModel* iomodel,int analysis_type);
void CreateNodes(Nodes*nodes, IoModel* iomodel,int analysis,int finite_element,int approximation=NoneApproximationEnum);

/*partitioning: */
void ElementsAndVerticesPartitioning(bool** pmy_elements, int** pmy_vertices, IoModel* iomodel);
void NodesPartitioning(bool** pmy_nodes,bool* my_elements, int* my_vertices,  IoModel* iomodel, bool continuous);
void FacesPartitioning(bool** pmy_faces,IoModel* iomodel);
void EdgesPartitioning(bool** pmy_edges,IoModel* iomodel);

/*Mesh properties*/
void CreateEdges(IoModel* iomodel);
void CreateFaces(IoModel* iomodel);
void CreateFaces3d(IoModel* iomodel);
void EdgeOnBoundaryFlags(bool** pflags,IoModel* iomodel);

/*Connectivity*/
void CreateSingleNodeToElementConnectivity(IoModel* iomodel);
void CreateNumberNodeToElementConnectivity(IoModel* iomodel);
#endif
