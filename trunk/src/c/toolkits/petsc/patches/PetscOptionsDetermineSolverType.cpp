/*!\file PetscOptionsDetermineSolverType.cpp: from the petsc options, determine what kind of solver
 * we are using.
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Petsc includes: */
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"

#include "./petscpatches.h"

#include "../../../shared/shared.h"

void PetscOptionsDetermineSolverType(int* psolver_type){

	char option[100];
   #if _PETSC_MAJOR_ < 3 || (_PETSC_MAJOR_ == 3 && _PETSC_MINOR_ < 2)
	PetscTruth flag;
	#else
	PetscBool flag;
	#endif

	/*output: */
	int solver_type=PETSCPACKAGE;

	/*retrieve mat_type option: */
	#if (_PETSC_MAJOR_==3) && (_PETSC_MINOR_>=7)
	PetscOptionsGetString(NULL,PETSC_NULL,"-mat_type",&option[0],100,&flag);
	#else
	PetscOptionsGetString(PETSC_NULL,"-mat_type",&option[0],100,&flag);
	#endif

	if (strcmp(option,"aijmumps")==0){
		solver_type=MUMPSPACKAGE_LU;
	}
	if (strcmp(option,"sbaijmumps")==0){
		solver_type=MUMPSPACKAGE_CHOL;
	}
	if (strcmp(option,"aijspooles")==0){
		solver_type=SPOOLESPACKAGE_LU;
	}
	if (strcmp(option,"sbaijspooles")==0){
		solver_type=SPOOLESPACKAGE_CHOL;
	}
	if (strcmp(option,"superlu_dist")==0){
		solver_type=SUPERLUDISTPACKAGE;
	}
	if (strcmp(option,"")==0){
		solver_type=SUPERLUDISTPACKAGE;
	}

	#if _PETSC_MAJOR_ >= 3 
		#if (_PETSC_MINOR_>=7)
		PetscOptionsGetString(NULL,PETSC_NULL,"-pc_factor_mat_solver_package",&option[0],100,&flag);
		#else
		PetscOptionsGetString(PETSC_NULL,"-pc_factor_mat_solver_package",&option[0],100,&flag);
		#endif
	if (strcmp(option,"mumps")==0){
		solver_type=MUMPSPACKAGE_LU;
	}
	#endif

	#if (_PETSC_MAJOR_==3) && (_PETSC_MINOR_>=7)
	PetscOptionsGetString(NULL,PETSC_NULL,"-issm_option_solver",&option[0],100,&flag);
	#else
	PetscOptionsGetString(PETSC_NULL,"-issm_option_solver",&option[0],100,&flag);
	#endif
	if (strcmp(option,"FS")==0){
		solver_type=FSSolverEnum;
	}

	*psolver_type=solver_type;
}
