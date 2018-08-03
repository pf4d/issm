/*!\file CreateJacobianMatrixx
 * \brief retrieve vector from inputs in elements
 */

#include "./CreateJacobianMatrixx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../AllocateSystemMatricesx/AllocateSystemMatricesx.h"

void CreateJacobianMatrixx(Matrix<IssmDouble>** pJff,FemModel* femmodel,IssmDouble kmax){

	int      i;
	int      configuration_type,analysisenum;
	Element *element = NULL;
	Load    *load    = NULL;
	Matrix<IssmDouble>* Jff = NULL;

	/*Checks*/
	_assert_(femmodel && femmodel->nodes && femmodel->elements);

	/*Recover some parameters*/
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	femmodel->parameters->FindParam(&analysisenum,AnalysisTypeEnum);
	Analysis* analysis = EnumToAnalysis(analysisenum);

	/*Initialize Jacobian Matrix*/
	AllocateSystemMatricesx(&Jff,NULL,NULL,NULL,femmodel);

	/*Create and assemble matrix*/
	for(i=0;i<femmodel->elements->Size();i++){
		element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ElementMatrix* Je = analysis->CreateJacobianMatrix(element);
		if(Je) Je->AddToGlobal(Jff);
		delete Je;
	}
	for (i=0;i<femmodel->loads->Size();i++){
		load=(Load*)femmodel->loads->GetObjectByOffset(i);
		if(load->InAnalysis(configuration_type)) load->CreateJacobianMatrix(Jff);
		if(load->InAnalysis(configuration_type)) load->PenaltyCreateJacobianMatrix(Jff,kmax);
	}
	Jff->Assemble();

	/*Assign output pointer*/
	delete analysis;
	*pJff=Jff;

}
