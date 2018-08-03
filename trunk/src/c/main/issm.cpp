/*!\file:  issm.cpp
 * \brief: ISSM main program
 */ 

#include "./issm.h"

int main(int argc,char **argv){

	/*Initialize exception trapping: */
	ExceptionTrapBegin();

	/*Initialize environment (MPI, PETSC, MUMPS, etc ...)*/
	ISSM_MPI_Comm comm_init=EnvironmentInit(argc,argv);

	/*Initialize femmodel from arguments provided command line: */
	FemModel *femmodel = new FemModel(argc,argv,comm_init);

	/*Solve: */
	femmodel->Solve();

	/*Output results: */
	OutputResultsx(femmodel);

	/*Wrap up: */
	femmodel->CleanUp();

	/*Delete Model: */
	delete femmodel;

	/*Finalize environment:*/
	EnvironmentFinalize();

	/*Finalize exception trapping: */
	ExceptionTrapEnd();

	/*Return unix success: */
	return 0; 
}
