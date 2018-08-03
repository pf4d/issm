/*!\file: solutionsequence_hydro_nonlinear.cpp
 * \brief: core of the hydro solution 
 */ 

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
/*FIXME, dirty hack to get the solutionsequence linear needed to compute the slopes*/
#include "../solutionsequences/solutionsequences.h"

void solutionsequence_hydro_nonlinear(FemModel* femmodel){
	/*solution : */
	Vector<IssmDouble>* ug_sed=NULL; 
	Vector<IssmDouble>* uf_sed=NULL; 
	Vector<IssmDouble>* uf_sed_sub_iter=NULL; 
	Vector<IssmDouble>* ug_sed_main_iter=NULL;

	Vector<IssmDouble>* ug_epl=NULL; 
	Vector<IssmDouble>* uf_epl=NULL;
	Vector<IssmDouble>* uf_epl_sub_iter=NULL; 
	Vector<IssmDouble>* ug_epl_sub_iter=NULL;
	Vector<IssmDouble>* ug_epl_main_iter=NULL;

	Vector<IssmDouble>* ys=NULL; 
	Vector<IssmDouble>* dug=NULL;
	Vector<IssmDouble>* duf=NULL;

	Matrix<IssmDouble>* Kff=NULL;
	Matrix<IssmDouble>* Kfs=NULL;

	Vector<IssmDouble>* pf=NULL;
	Vector<IssmDouble>* df=NULL;

	HydrologyDCInefficientAnalysis* inefanalysis = NULL;
	HydrologyDCEfficientAnalysis* effanalysis = NULL;
	
	bool       sedconverged,eplconverged,hydroconverged;
	bool       isefficientlayer;
	int        constraints_converged;
	int        num_unstable_constraints;
	int        sedcount,eplcount,hydrocount;
	int        hydro_maxiter;
	IssmDouble sediment_kmax;
	IssmDouble eps_hyd;
	IssmDouble ndu_sed,nu_sed;
	IssmDouble ndu_epl,nu_epl;
	IssmDouble ThickCount,L2Count;
	
	/*Recover parameters: */
	femmodel->SetCurrentConfiguration(HydrologyDCInefficientAnalysisEnum);
	femmodel->parameters->FindParam(&isefficientlayer,HydrologydcIsefficientlayerEnum);
	femmodel->parameters->FindParam(&hydro_maxiter,HydrologydcMaxIterEnum);
	femmodel->parameters->FindParam(&eps_hyd,HydrologydcRelTolEnum);

	hydrocount=1;
	hydroconverged=false;
	/*We don't need the outer loop if only one layer is used*/
	if(!isefficientlayer) hydroconverged=true;

	/*Retrieve inputs as the initial state for the non linear iteration*/
	GetSolutionFromInputsx(&ug_sed,femmodel);	
	Reducevectorgtofx(&uf_sed, ug_sed, femmodel->nodes,femmodel->parameters);

	if(isefficientlayer) {
		inefanalysis = new HydrologyDCInefficientAnalysis();
		effanalysis = new HydrologyDCEfficientAnalysis();
		femmodel->SetCurrentConfiguration(HydrologyDCEfficientAnalysisEnum);
		GetSolutionFromInputsx(&ug_epl,femmodel);
		/*Initialize the element mask*/
		inefanalysis->ElementizeEplMask(femmodel);
		effanalysis->InitZigZagCounter(femmodel);
	}
	/*The real computation starts here, outermost loop is on the two layer system*/
	for(;;){

		sedcount=1;
		eplcount=1;

		/*If there is two layers we need an outer loop value to compute convergence*/
		if(isefficientlayer){
			ug_sed_main_iter=ug_sed->Duplicate();
			ug_sed->Copy(ug_sed_main_iter);
			ug_epl_main_iter=ug_epl->Duplicate();
			ug_epl->Copy(ug_epl_main_iter);
		}
		/*Loop on sediment layer to deal with transfer and head value*/
		femmodel->SetCurrentConfiguration(HydrologyDCInefficientAnalysisEnum);
		InputUpdateFromConstantx(femmodel,true,ResetPenaltiesEnum);
		InputUpdateFromConstantx(femmodel,false,ConvergedEnum);
		femmodel->UpdateConstraintsx();
		femmodel->parameters->SetParam(HydrologySedimentEnum,HydrologyLayerEnum);
		
		/*Reset constraint on the ZigZag Lock*/
		ResetConstraintsx(femmodel);

		/*{{{*//*Treating the sediment*/
		for(;;){
			sedconverged=false;
			uf_sed_sub_iter=uf_sed->Duplicate();_assert_(uf_sed_sub_iter);
			uf_sed->Copy(uf_sed_sub_iter);
			/*{{{*//*Loop on the sediment layer to deal with the penalization*/
			for(;;){
				/*{{{*/ /*Core of the computation*/
				if(VerboseSolution()) _printf0_("Building Sediment Matrix...\n");
				SystemMatricesx(&Kff,&Kfs,&pf,&df,&sediment_kmax,femmodel);
				CreateNodalConstraintsx(&ys,femmodel->nodes,HydrologyDCInefficientAnalysisEnum);
				Reduceloadx(pf,Kfs,ys); delete Kfs;
				delete uf_sed;
				Solverx(&uf_sed,Kff,pf,uf_sed_sub_iter,df,femmodel->parameters);
				delete Kff; delete pf; delete df;
				delete ug_sed;
				Mergesolutionfromftogx(&ug_sed,uf_sed,ys,femmodel->nodes,femmodel->parameters); delete ys;
				InputUpdateFromSolutionx(femmodel,ug_sed);
				ConstraintsStatex(&constraints_converged,&num_unstable_constraints,femmodel);
				/*}}}*/
				if (!sedconverged){
					if(VerboseConvergence()) _printf0_("   # Sediment unstable constraints = " << num_unstable_constraints << "\n");
					if(num_unstable_constraints==0) sedconverged = true;
					if (sedcount>=hydro_maxiter){
						_error_("   maximum number of Sediment iterations (" << hydro_maxiter << ") exceeded");
					}
				}
				/*Add an iteration and get out of the loop if the penalisation is converged*/
				sedcount++;
				if(sedconverged)break;
			}
		
			/*}}}*//*End of the sediment penalization loop*/
			sedconverged=false;
			
			/*Checking convegence on the value of the sediment head*/
			duf=uf_sed_sub_iter->Duplicate();_assert_(duf);
			uf_sed_sub_iter->Copy(duf);
			duf->AYPX(uf_sed,-1.0);
			ndu_sed=duf->Norm(NORM_TWO);
			delete duf;
			nu_sed=uf_sed_sub_iter->Norm(NORM_TWO);
			if (xIsNan<IssmDouble>(ndu_sed) || xIsNan<IssmDouble>(nu_sed)) _error_("convergence criterion is NaN!");
			if (ndu_sed==0.0 && nu_sed==0.0) nu_sed=1.0e-6; /*Hacking the case where the layer is empty*/
			if(VerboseConvergence()) _printf0_(setw(50) << left << "   Inner Sediment Convergence criterion:" << ndu_sed/nu_sed*100 << "%, aiming lower than " << eps_hyd*10*100 << " %\n");
			if((ndu_sed/nu_sed)<eps_hyd*10.){
				if(VerboseConvergence()) _printf0_("   # Inner sediment convergence achieve \n");
				sedconverged=true;
			}
			delete uf_sed_sub_iter;
			if(sedconverged){
				femmodel->parameters->SetParam(sediment_kmax,HydrologySedimentKmaxEnum);
				InputUpdateFromConstantx(femmodel,sedconverged,ConvergedEnum);
				InputUpdateFromSolutionx(femmodel,ug_sed);
				InputUpdateFromConstantx(femmodel,sediment_kmax,HydrologySedimentKmaxEnum);
				break;
			}
		}
		/*}}}*//*End of the global sediment loop*/
		/*{{{*//*Now dealing with the EPL in the same way*/
		if(isefficientlayer){
			femmodel->SetCurrentConfiguration(HydrologyDCEfficientAnalysisEnum);
			/*updating mask*/
			if(VerboseSolution()) _printf0_("==updating mask...\n");
			femmodel->HydrologyEPLupdateDomainx(&ThickCount);
			inefanalysis->ElementizeEplMask(femmodel);
			InputUpdateFromConstantx(femmodel,true,ResetPenaltiesEnum);
			InputUpdateFromConstantx(femmodel,false,ConvergedEnum);
			femmodel->parameters->SetParam(HydrologyEfficientEnum,HydrologyLayerEnum);

			for(;;){
				eplconverged=false;
				ug_epl_sub_iter=ug_epl->Duplicate();_assert_(ug_epl_sub_iter);
				ug_epl->Copy(ug_epl_sub_iter);
				/*{{{*//*Retrieve the EPL head slopes and compute EPL Thickness*/
				if(VerboseSolution()) _printf0_("computing EPL Head slope...\n");
				femmodel->SetCurrentConfiguration(L2ProjectionEPLAnalysisEnum);
				femmodel->UpdateConstraintsL2ProjectionEPLx(&L2Count);
				inefanalysis->ElementizeEplMask(femmodel);
				femmodel->parameters->SetParam(EplHeadSlopeXEnum,InputToL2ProjectEnum);
				solutionsequence_linear(femmodel);
				femmodel->parameters->SetParam(EplHeadSlopeYEnum,InputToL2ProjectEnum);
				solutionsequence_linear(femmodel);

				femmodel->SetCurrentConfiguration(HydrologyDCEfficientAnalysisEnum);
				effanalysis->ComputeEPLThickness(femmodel);
				//updating mask after the computation of the epl thickness (Allow to close too thin EPL)
				femmodel->HydrologyEPLupdateDomainx(&ThickCount);
				inefanalysis->ElementizeEplMask(femmodel);
				/*}}}*/
					
				if(VerboseSolution()) _printf0_("Building EPL Matrix...\n");
				SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
				CreateNodalConstraintsx(&ys,femmodel->nodes,HydrologyDCEfficientAnalysisEnum);
				Reduceloadx(pf,Kfs,ys); delete Kfs;
				delete uf_epl;
				Solverx(&uf_epl,Kff,pf,uf_epl_sub_iter,df,femmodel->parameters);
				delete Kff; delete pf; delete df;
				delete uf_epl_sub_iter;
				uf_epl_sub_iter=uf_epl->Duplicate();
				uf_epl->Copy(uf_epl_sub_iter);
				delete ug_epl; 
				Mergesolutionfromftogx(&ug_epl,uf_epl,ys,femmodel->nodes,femmodel->parameters); delete ys;
				InputUpdateFromSolutionx(femmodel,ug_epl);
				ConstraintsStatex(&constraints_converged,&num_unstable_constraints,femmodel);
						
				dug=ug_epl_sub_iter->Duplicate();_assert_(dug);
				ug_epl_sub_iter->Copy(dug);
				dug->AYPX(ug_epl,-1.0);
				ndu_epl=dug->Norm(NORM_TWO);
				delete dug;
				nu_epl=ug_epl_sub_iter->Norm(NORM_TWO);
				if (xIsNan<IssmDouble>(ndu_epl) || xIsNan<IssmDouble>(nu_epl)) _error_("convergence criterion is NaN!");
				if (ndu_epl==0.0 && nu_epl==0.0) nu_epl=1.0e-6; /*Hacking the case where the EPL is used but empty*/
				if(VerboseConvergence()) _printf0_(setw(50) << left << "   Inner EPL Convergence criterion:" << ndu_epl/nu_epl*100 << "%, aiming lower than " << eps_hyd*10*100 << " %\n");
				if((ndu_epl/nu_epl)<eps_hyd*10.) eplconverged=true;
				if (eplcount>=hydro_maxiter){
					_error_("   maximum number of EPL iterations (" << hydro_maxiter << ") exceeded");
				}
				//If there is some colapse go through sediment again
				/* if(ThickCount<L2Count)eplconverged=true; */
				eplcount++;
				
				delete ug_epl_sub_iter;
				if(eplconverged){
					if(VerboseSolution()) _printf0_("eplconverged...\n");
					InputUpdateFromConstantx(femmodel,eplconverged,ConvergedEnum);
					InputUpdateFromSolutionx(femmodel,ug_epl);
					effanalysis->ResetCounter(femmodel);
					break;
				}
			}
		}
		/*}}}*/ /*End of the global EPL loop*/

		/*{{{*/ /*Now dealing with the convergence of the whole system*/
		if(!hydroconverged){
			//compute norm(du)/norm(u)
			dug=ug_sed_main_iter->Duplicate(); _assert_(dug);
			ug_sed_main_iter->Copy(dug);	
			dug->AYPX(ug_sed,-1.0);
			ndu_sed=dug->Norm(NORM_TWO); 
			delete dug;
			nu_sed=ug_sed_main_iter->Norm(NORM_TWO);
			delete ug_sed_main_iter;
			if (xIsNan<IssmDouble>(ndu_sed) || xIsNan<IssmDouble>(nu_sed)) _error_("Sed convergence criterion is NaN!");
			if (ndu_sed==0.0 && nu_sed==0.0) nu_sed=1.0e-6; /*Hacking the case where the Sediment is used but empty*/
			dug=ug_epl_main_iter->Duplicate();_assert_(dug); 
			ug_epl_main_iter->Copy(dug); 
			dug->AYPX(ug_epl,-1.0);
			ndu_epl=dug->Norm(NORM_TWO); 
			delete dug;
			nu_epl=ug_epl_main_iter->Norm(NORM_TWO);
			delete ug_epl_main_iter;
			if (xIsNan<IssmDouble>(ndu_epl) || xIsNan<IssmDouble>(nu_epl)) _error_("EPL convergence criterion is NaN!");
			if (ndu_epl==0.0 && nu_epl==0.0) nu_epl=1.0e-6; /*Hacking the case where the EPL is used but empty*/
			if (!xIsNan<IssmDouble>(eps_hyd)){
				if ((ndu_epl/nu_epl)<eps_hyd && (ndu_sed/nu_sed)<(eps_hyd)){
					if (VerboseConvergence()) _printf0_(setw(50) << left << "   Converged after, " << hydrocount << " iterations \n");
					hydroconverged=true;
				}
				else{ 
					if(VerboseConvergence()) _printf0_(setw(50) << left << "   Sediment Convergence criterion:" << ndu_sed/nu_sed*100 << "%, aiming lower than " << eps_hyd*100 << " %\n");
					if(VerboseConvergence()) _printf0_(setw(50) << left << "   EPL Convergence criterion:" << ndu_epl/nu_epl*100 << "%, aiming lower than " << eps_hyd*100 << " %\n");
					hydroconverged=false;
				}
			}
			else _printf0_(setw(50) << left << "   Convergence criterion:" << ndu_sed/nu_sed*100 << " %\n");
			if (hydrocount>=hydro_maxiter){
				_error_("   maximum number for hydrological global iterations (" << hydro_maxiter << ") exceeded");
			}
		}
		hydrocount++;
		if(hydroconverged)break;
	}
	/*}}}*/
	if(isefficientlayer)InputUpdateFromSolutionx(femmodel,ug_epl);
	femmodel->SetCurrentConfiguration(HydrologyDCInefficientAnalysisEnum);
	InputUpdateFromSolutionx(femmodel,ug_sed);
	/*Free ressources: */
	delete ug_epl;
	delete ug_sed;
	delete uf_sed;
	delete uf_epl;
	delete uf_epl_sub_iter;
	delete inefanalysis;
	delete effanalysis;
}
