/*!\file:  KSPFree.cpp
 * \brief wrapper to KSPDestroy
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Petsc includes: */
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>

void KSPFree(KSP* pksp){

	#if _PETSC_MAJOR_ < 3 || (_PETSC_MAJOR_ == 3 && _PETSC_MINOR_ < 2)
	if(*pksp)KSPDestroy(*pksp);
	*pksp=NULL;
	#else
	if(*pksp)KSPDestroy(pksp);
	*pksp=NULL;
	#endif

}
