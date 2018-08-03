/*!\file: solutionsequence_linear.cpp
 * \brief: numerical core of linear solutions
 */ 

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_linear(FemModel* femmodel){

	/*intermediary: */
	Matrix<IssmDouble>*  Kff = NULL;
	Matrix<IssmDouble>*  Kfs = NULL;
	Vector<IssmDouble>*  ug  = NULL;
	Vector<IssmDouble>*  uf  = NULL;
	Vector<IssmDouble>*  pf  = NULL;
	Vector<IssmDouble>*  df  = NULL;
	Vector<IssmDouble>*  ys  = NULL;
	int  configuration_type;
	IssmDouble solver_residue_threshold;
	
	/*solver convergence test: */
	IssmDouble nKUF;
	IssmDouble nF;
	IssmDouble solver_residue;
	Vector<IssmDouble>* KU=NULL;
	Vector<IssmDouble>* KUF=NULL;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	femmodel->parameters->FindParam(&solver_residue_threshold,SettingsSolverResidueThresholdEnum);
	femmodel->UpdateConstraintsx();
	SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
	CreateNodalConstraintsx(&ys,femmodel->nodes,configuration_type);
	Reduceloadx(pf, Kfs, ys); delete Kfs;
	Solverx(&uf, Kff, pf, NULL, df, femmodel->parameters); 
	
	/*Check that the solver converged nicely: */
		
	//compute KUF = KU - F = K*U - F
	KU=uf->Duplicate(); Kff->MatMult(uf,KU);
	KUF=KU->Duplicate(); KU->Copy(KUF); KUF->AYPX(pf,-1.0);

	//compute norm(KUF), norm(F) and residue
	nKUF=KUF->Norm(NORM_TWO);
	nF=pf->Norm(NORM_TWO);
	solver_residue=nKUF/nF;

	if(!xIsNan(solver_residue_threshold) && solver_residue>solver_residue_threshold)_error_("   solver residue too high!: norm(KU-F)/norm(F)=" << solver_residue << "\n");


	//clean up
	delete KU; delete KUF;
	delete Kff; delete pf; delete df;

	//#ifdef  _HAVE_ADOLC_
	//        for (int i =0; i<uf->svector->M; ++i) {
	//          ADOLC_DUMP_MACRO(uf->svector->vector[i]);
	//        }
	//#endif
	
	Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters);delete uf; delete ys;
	InputUpdateFromSolutionx(femmodel,ug); 
	delete ug;  
}
