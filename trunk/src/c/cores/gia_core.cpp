/*!\file: gia_core.cpp
 * \brief: core of the GIA solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"
void gia_core(FemModel* femmodel){

	Vector<IssmDouble> *wg    = NULL;
	Vector<IssmDouble> *dwdtg = NULL;
	IssmDouble          *x    = NULL;
	IssmDouble          *y    = NULL;

	/*parameters: */
	bool save_results;
	int  gsize;
	int  configuration_type;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);

	if(VerboseSolution()) _printf0_("   computing GIA\n");

	/*Call on core computations: */
	femmodel->SetCurrentConfiguration(GiaIvinsAnalysisEnum);

	/*Figure out size of g-set deflection vector and allocate solution vector: */
	gsize      = femmodel->nodes->NumberOfDofs(configuration_type,GsetEnum);
	wg = new Vector<IssmDouble>(gsize);
	dwdtg = new Vector<IssmDouble>(gsize);

	/*first, recover x and y vectors from vertices: */
	VertexCoordinatesx(&x,&y,NULL,femmodel->vertices); //no need for z coordinate

	/*call the main module: */
	femmodel->Deflection(wg,dwdtg,x,y);

	/*assemble vector: */
	wg->Assemble();
	dwdtg->Assemble();

	InputUpdateFromVectorx(femmodel,wg,GiaWEnum,VertexSIdEnum);
	InputUpdateFromVectorx(femmodel,dwdtg,GiadWdtEnum,VertexSIdEnum);

	if(save_results){
		if(VerboseSolution()) _printf0_("   saving results\n");
		int outputs[2] = {GiaWEnum,GiadWdtEnum};
		femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],2);
	}
	
	xDelete<IssmDouble>(x);
	xDelete<IssmDouble>(y);


}
