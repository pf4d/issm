/*!\file ModelProcessorx
 * \brief: create datasets using input binary file and a set of requested analyses
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void ModelProcessorx(Elements** pelements, Nodes** pnodes, Vertices** pvertices, Materials** pmaterials, Constraints** pconstraints, Loads** ploads, Parameters** pparameters,IoModel* iomodel,FILE* toolkitfile, char* rootpath,const int solution_enum,const int nummodels,const int* analysis_enum_list){

	int   i,analysis_enum,verbose;

	/*Initialize datasets*/
	Elements    *elements    = new Elements();
	Nodes       *nodes       = new Nodes();
	Vertices    *vertices    = new Vertices();
	Materials   *materials   = new Materials();
	Constraints *constraints = new Constraints();
	Loads       *loads       = new Loads();
	Parameters  *parameters  = new Parameters();

	/*Fetch parameters: */
	iomodel->FindConstant(&verbose,"md.verbose");
	SetVerbosityLevel(verbose);

	if(VerboseMProcessor()) _printf0_("   starting model processor \n");

	/*Partition Elements and Nodes*/
	ElementsAndVerticesPartitioning(&iomodel->my_elements,&iomodel->my_vertices,iomodel);

	/*Create elements, vertices and materials, independent of analysis_enum: */
	CreateElementsVerticesAndMaterials(elements,vertices,materials,iomodel,nummodels);

	/*Create Parameters*/
	CreateParameters(parameters,iomodel,rootpath,toolkitfile,solution_enum);

	for(i=0;i<nummodels;i++){

		analysis_enum=analysis_enum_list[i];
		parameters->AddObject(new IntParam(AnalysisCounterEnum,i));

		if(VerboseMProcessor()) _printf0_("   creating datasets for analysis " << EnumToStringx(analysis_enum) << "\n");
		Analysis* analysis = EnumToAnalysis(analysis_enum);
		analysis->UpdateParameters(parameters,iomodel,solution_enum,analysis_enum);
		analysis->CreateNodes(nodes,iomodel);
		analysis->CreateConstraints(constraints,iomodel);
		analysis->CreateLoads(loads,iomodel);
		analysis->UpdateElements(elements,iomodel,i,analysis_enum);
		delete analysis;


		/* Update counters, because we have created more nodes, loads and
		 * constraints, and ids for objects created in next call to CreateDataSets
		 * will need to start at the end of the updated counters: */
		if(nodes->Size()) iomodel->nodecounter = nodes->MaximumId();
		iomodel->loadcounter       = loads->NumberOfLoads();
		iomodel->constraintcounter = constraints->NumberOfConstraints();

		/*Make sure nodecounter is at least 0 (if no node exists, maxid will be -1*/
		_assert_(iomodel->nodecounter>=0);
	}

	/*Solution specific updates*/
	UpdateElementsAndMaterialsControl(elements,materials,iomodel);
	#ifdef _HAVE_DAKOTA_
	UpdateElementsAndMaterialsDakota(elements,materials,iomodel);
	#endif
	if(solution_enum==TransientSolutionEnum){
		UpdateElementsTransient(elements,parameters,iomodel,analysis_enum);
	}

	/*Output definitions dataset: */
	CreateOutputDefinitions(elements,parameters,iomodel);

	/* Sort datasets:
	 * All our datasets are already ordered by ids. Set presort flag so that
	 * later on, when sorting is requested on these datasets, it will not be
	 * redone: */
	elements->Presort();
	nodes->Presort();
	vertices->Presort();
	loads->Presort();
	materials->Presort();

	constraints->Presort();
	if(VerboseMProcessor()) _printf0_("   done with model processor \n");

	/*Assign output pointers:*/
	*pelements    = elements;
	*pnodes       = nodes;
	*pvertices    = vertices;
	*pmaterials   = materials;
	*pconstraints = constraints;
	*ploads       = loads;
	*pparameters  = parameters;
}
