/*!\file SolverxPetsc
 * \brief Petsc implementation of solver
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./PetscSolver.h"
#include "../../../shared/Numerics/Verbosity.h"
#include "../../../shared/MemOps/MemOps.h"
#include "../../../shared/Exceptions/exceptions.h"
#include "../../../shared/io/Comm/IssmComm.h"
#include "../../../shared/Enum/Enum.h"

void	PetscSolve(PetscVec** puf, PetscMat* Kff, PetscVec* pf, PetscVec* uf0,PetscVec* df, Parameters* parameters){ /*{{{*/

	PetscVec* uf=new PetscVec();

	Vec uf0_vector = NULL;
	Vec df_vector  = NULL;

	if(uf0) uf0_vector = uf0->vector;
	if(df)  df_vector  = df->vector;

	SolverxPetsc(&uf->vector, Kff->matrix, pf->vector, uf0_vector, df_vector, parameters);

	*puf=uf;

}
/*}}}*/
void	SolverxPetsc(Vec* puf, Mat Kff, Vec pf, Vec uf0,Vec df, Parameters* parameters){ /*{{{*/

	/*Output: */
	Vec        uf               = NULL;

	/*Intermediary: */
	int        local_m,local_n,global_m,global_n;

	/*Solver */
	KSP        ksp              = NULL;
	PC         pc               = NULL;
	int        iteration_number;
	int        solver_type;
	bool       fromlocalsize    = true;
	#if _PETSC_MAJOR_ < 3 || (_PETSC_MAJOR_ == 3 && _PETSC_MINOR_ < 2)
	PetscTruth flag,flg;
	#else
	PetscBool flag,flg;
	#endif

	/*FS: */
	IS         isv=NULL;
	IS         isp=NULL;

	#if _PETSC_MAJOR_ >= 3 
	char ksp_type[50];
	#endif

	/*Display message*/
	#if _PETSC_MAJOR_ < 3 || (_PETSC_MAJOR_ == 3 && _PETSC_MINOR_ < 2)
	if(VerboseSolver())PetscOptionsPrint(stdout);
	#else
		#if _PETSC_MINOR_<7
		if(VerboseSolver())PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD);
		#else
		if(VerboseSolver())PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD);
		#endif
	#endif

	/*First, check that f-set is not NULL, i.e. model is fully constrained:*/ 
	_assert_(Kff);
	MatGetSize(Kff,&global_m,&global_n); _assert_(global_m==global_n);
	if(!global_n){
		*puf=NewVec(0,IssmComm::GetComm()); return;
	}

	/*Initial guess */
	/*Now, check that we are not giving an initial guess to the solver, if we are running a direct solver: */
	#if _PETSC_MAJOR_ >= 3 
		#if _PETSC_MINOR_<7
		PetscOptionsGetString(PETSC_NULL,"-ksp_type",ksp_type,49,&flg);
		#else
		PetscOptionsGetString(NULL,PETSC_NULL,"-ksp_type",ksp_type,49,&flg);
		#endif
	if (strcmp(ksp_type,"preonly")==0)uf0=NULL;
	#endif

	/*If initial guess for the solution exists, use it to create uf, otherwise, 
	 * duplicate the right hand side so that the solution vector has the same structure*/
	if(uf0){
		VecDuplicate(uf0,&uf); VecCopy(uf0,uf);
	}
	else{
		MatGetLocalSize(Kff,&local_m,&local_n);uf=NewVec(local_n,IssmComm::GetComm(),fromlocalsize);
	}

	/*Process petsc options to see if we are using special types of external solvers*/
	PetscOptionsDetermineSolverType(&solver_type);

	/*Check the solver is available*/
	if(solver_type==MUMPSPACKAGE_LU || solver_type==MUMPSPACKAGE_CHOL){
		#if _PETSC_MAJOR_ >=3
			#ifndef _HAVE_MUMPS_
			_error_("requested MUMPS solver, which was not compiled into ISSM!\n");
			#endif
		#endif
	}

	/*Prepare solver*/
	KSPCreate(IssmComm::GetComm(),&ksp);
	#if (_PETSC_MAJOR_==3) && (_PETSC_MINOR_>=5)
		KSPSetOperators(ksp,Kff,Kff);
	#else
		KSPSetOperators(ksp,Kff,Kff,DIFFERENT_NONZERO_PATTERN);
	#endif
	KSPSetFromOptions(ksp);

	#if _PETSC_MAJOR_==3
	/*Specific solver?: */
	KSPGetPC(ksp,&pc);
	if (solver_type==MUMPSPACKAGE_LU){
		#if _PETSC_MINOR_==1
		PCFactorSetMatSolverPackage(pc,MAT_SOLVER_MUMPS);
		#else
		PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
		#endif
	}

	/*FS: */
	if (solver_type==FSSolverEnum){
		/*Make indices out of doftypes: */
		if(!df)_error_("need doftypes for FS solver!\n");
		DofTypesToIndexSet(&isv,&isp,df,FSSolverEnum);

		/*Set field splits: */
		KSPGetPC(ksp,&pc);
		#if _PETSC_MINOR_==1
		PCFieldSplitSetIS(pc,isv);
		PCFieldSplitSetIS(pc,isp);
		#else
		PCFieldSplitSetIS(pc,PETSC_NULL,isv);
		PCFieldSplitSetIS(pc,PETSC_NULL,isp);
		#endif

	}
	#endif

	/*If there is an initial guess for the solution, use it
	 * except if we are using the MUMPS direct solver
	 * where any initial guess will crash Petsc*/
	if (uf0){
		if((solver_type!=MUMPSPACKAGE_LU) && (solver_type!=MUMPSPACKAGE_CHOL) && (solver_type!=SPOOLESPACKAGE_LU) && (solver_type!=SPOOLESPACKAGE_CHOL) && (solver_type!=SUPERLUDISTPACKAGE)){
			KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
		}
	}

	/*Solve: */
	if(VerboseSolver())KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
	KSPSolve(ksp,pf,uf);

	/*Check convergence*/
	KSPGetIterationNumber(ksp,&iteration_number);
	if (iteration_number<0) _error_("Solver diverged at iteration number: " << -iteration_number);

	/*Free resources:*/
	KSPFree(&ksp);

	/*Assign output pointers:*/
	*puf=uf;
} 
/*}}}*/
void DofTypesToIndexSet(IS* pisv, IS* pisp, Vec df,int typeenum){ /*{{{*/

	/*output: */
	IS          isv=NULL;
	IS          isp=NULL;

	int         start,end;
	IssmDouble*     df_local=NULL;
	int         df_local_size;

	int*     pressure_indices=NULL;
	int*     velocity_indices=NULL;
	int      pressure_num=0;
	int      velocity_num=0;
	int      pressure_count=0;
	int      velocity_count=0;

	if(typeenum==FSSolverEnum){

		/*Ok, recover doftypes vector values and indices: */
		VecGetOwnershipRange(df,&start,&end);
		VecGetLocalSize(df,&df_local_size);
		VecGetArray(df,&df_local);

		pressure_num=0;
		velocity_num=0;
		for(int i=0;i<df_local_size;i++){
			if (df_local[i]==PressureEnum)pressure_num++;
			else velocity_num++;
		}

		/*Allocate indices: */
		if(pressure_num)pressure_indices=xNew<int>(pressure_num);
		if(velocity_num)velocity_indices=xNew<int>(velocity_num);

		pressure_count=0;
		velocity_count=0;
		for(int i=0;i<df_local_size;i++){
			if (df_local[i]==PressureEnum){
				pressure_indices[pressure_count]=start+i;
				pressure_count++;
			}
			if (df_local[i]==VelocityEnum){
				velocity_indices[velocity_count]=start+i;
				velocity_count++;
			}
		}
		VecRestoreArray(df,&df_local);

		/*Create indices sets: */
		#if _PETSC_MAJOR_ < 3 || (_PETSC_MAJOR_ == 3 && _PETSC_MINOR_ < 2)
		ISCreateGeneral(IssmComm::GetComm(),pressure_num,pressure_indices,&isp);
		ISCreateGeneral(IssmComm::GetComm(),velocity_num,velocity_indices,&isv);
		#else
		ISCreateGeneral(IssmComm::GetComm(),pressure_num,pressure_indices,PETSC_COPY_VALUES,&isp);
		ISCreateGeneral(IssmComm::GetComm(),velocity_num,velocity_indices,PETSC_COPY_VALUES,&isv);
		#endif
	}

	/*Free ressources:*/
	xDelete<int>(pressure_indices);
	xDelete<int>(velocity_indices);

	/*Assign output pointers:*/
	*pisv=isv;
	*pisp=isp;
}
/*}}}*/
