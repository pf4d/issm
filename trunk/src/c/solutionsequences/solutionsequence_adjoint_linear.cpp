/*!\file: solutionsequence_adjoint_linear.cpp
 * \brief: numerical core of linear solutions
 */ 

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_adjoint_linear(FemModel* femmodel){
	/*This is axactly the same solutionsequence as solutionsequence_linear except that Reduceloadfromgtofx and Mergesolutionfromftogx
	 * use the flag "true" so that all spc are taken as 0*/

	/*intermediary: */
	Matrix<IssmDouble>*  Kff = NULL;
	Matrix<IssmDouble>*  Kfs = NULL;
	Vector<IssmDouble>*  ug  = NULL;
	Vector<IssmDouble>*  uf  = NULL;
	Vector<IssmDouble>*  pf  = NULL;
	Vector<IssmDouble>*  df  = NULL;
	Vector<IssmDouble>*  ys  = NULL;
	int  configuration_type;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	femmodel->UpdateConstraintsx();

	SystemMatricesx(&Kff, &Kfs, &pf, &df, NULL,femmodel);
	CreateNodalConstraintsx(&ys,femmodel->nodes,configuration_type);
	Reduceloadx(pf, Kfs, ys,true); delete Kfs; //true means spc = 0
	Solverx(&uf, Kff, pf, NULL, df, femmodel->parameters); delete Kff; delete pf; delete df;
	Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters,true); delete ys; //true means spc0
	InputUpdateFromSolutionx(femmodel,ug);
	delete ug; delete uf;
}
