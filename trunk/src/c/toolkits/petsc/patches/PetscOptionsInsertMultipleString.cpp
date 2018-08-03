/*!\file PetscOptionsInsertMultipleString.cpp
 * \brief: create distributed Petsc vector.
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

#include "./petscpatches.h"

#include "../../../shared/shared.h"

void PetscOptionsInsertMultipleString(char* options_string){

	/*The list of options is going to be pairs of the type "-option option_value"*/
	#if _PETSC_MAJOR_ == 2
		PetscToken *token=NULL ;
	#else
		PetscToken token=NULL;
	#endif
	char* first=NULL;
	char* second=NULL;
	size_t len;
	int first_token=1;

	PetscTokenCreate(options_string,' ',&token);
	for (;;){

		/*Read next tokens*/
		if(first_token){
			PetscTokenFind(token,&first);
		}
		PetscTokenFind(token,&second);

		if (!first){
			/*We are at the end of options*/
			break;
		}
		if(!second){
			/*We have no second value, just take first
			 * and set the option, then end the token analysis.*/
			if(first[0]!='-'){
				/*This is not good, the option does not have '-'! Get out*/
				_error_("Option " << first << " should be preceded by '-'!");
			}
			/*Reduce first to bare option value*/
			PetscStrlen(first,&len);
			while (len > 0 && first[len-1] == ' ') {
				len--; first[len] = 0;
			}
			#if _PETSC_MAJOR_ == 3 && _PETSC_MINOR_ < 7
			PetscOptionsSetValue(first,second);
			#else
			PetscOptionsSetValue(NULL,first,second);
			#endif
			break;
		}
		else{
			/*Ok, we have a second token coming in. Is it another option, or 'first' option's value?*/
			if (second[0]=='-'){
				/*Second is another option, ignore it*/
				PetscStrlen(first,&len);
				while (len > 0 && first[len-1] == ' ' ) {
					len--; first[len] = 0;
				}
				#if _PETSC_MAJOR_ == 3 && _PETSC_MINOR_ < 7
				PetscOptionsSetValue(first,NULL);
				#else
				PetscOptionsSetValue(NULL,first,NULL);
				#endif
				/*Preparing next loop step*/
				first=second;
				first_token=0;
			}
			else{
				/*Second is 'first' option's value*/
				PetscStrlen(second,&len);
				while (len > 0 && second[len-1] == ' ') {
					len--; second[len] = 0;
				}
				#if _PETSC_MAJOR_ == 3 && _PETSC_MINOR_ < 7
				PetscOptionsSetValue(first,second);
				#else
				PetscOptionsSetValue(NULL,first,second);
				#endif
				first_token=1;
			}
		}
	}

#if _PETSC_MAJOR_ == 3 && _PETSC_MINOR_ > 2
	PetscTokenDestroy(&token);
#else
	PetscTokenDestroy(token);
#endif
}
