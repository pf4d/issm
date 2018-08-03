/*!\file Solverx
 * \brief solver
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Solverx.h"
#include "../../shared/shared.h"

void	Solverx(Vector<IssmDouble>** puf, Matrix<IssmDouble>* Kff, Vector<IssmDouble>* pf, Vector<IssmDouble>* uf0,Vector<IssmDouble>* df, Parameters* parameters){

	/*intermediary: */
	Solver<IssmDouble> *solver=NULL;

	/*output: */
	Vector<IssmDouble> *uf=NULL;

	if(VerboseModule()) _printf0_("   Solving matrix system\n");

	/*Initialize solver: */
	solver=new Solver<IssmDouble>(Kff,pf,uf0,df,parameters);

	/*Solve:*/
	uf=solver->Solve();

	/*clean up and assign output pointers:*/
	delete solver;
	*puf=uf;
}
