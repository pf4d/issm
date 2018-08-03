/*!\file GetSolutionFromInputsx
 * \brief: update datasets using  parameter inputs
 */

#include "./GetSolutionFromInputsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void GetSolutionFromInputsx(Vector<IssmDouble>** psolution,FemModel* femmodel){

	/*intermediary: */
	int      gsize;
	int      configuration,analysisenum;

	/*output: */
	Vector<IssmDouble>* solution=NULL;

	if(VerboseModule()) _printf0_("   Get solution from inputs\n");

	/*retrive parameters: */
	femmodel->parameters->FindParam(&configuration,ConfigurationTypeEnum);
	femmodel->parameters->FindParam(&analysisenum,AnalysisTypeEnum);

	/*Get size of vector: */
	gsize=femmodel->nodes->NumberOfDofs(configuration,GsetEnum);
	if(gsize==0) _error_("Allocating a Vec of size 0 as gsize=0 for configuration "<<EnumToStringx(configuration));

	/*Initialize solution: */
	solution=new Vector<IssmDouble>(gsize);

	/*Go through elements and plug solution: */
	Analysis* analysis = EnumToAnalysis(analysisenum);
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		analysis->GetSolutionFromInputs(solution,element);
	}
	delete analysis;

	/*Assemble vector: */
	solution->Assemble();

	/*Assign output pointers:*/
	*psolution=solution;
}
