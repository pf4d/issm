/*!\file:  MatFree.cpp
 * \brief wrapper to MatDestroy
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Petsc includes: */
#include <petscmat.h>
#include <petscmat.h>
#include <petscksp.h>

void MatFree(Mat* pmat){

	#if _PETSC_MAJOR_ < 3 || (_PETSC_MAJOR_ == 3 && _PETSC_MINOR_ < 2)
	if(*pmat)MatDestroy(*pmat);
	*pmat=NULL;
	#else
	if(*pmat)MatDestroy(pmat);
	*pmat=NULL;
	#endif

}
