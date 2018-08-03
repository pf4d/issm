/*
 * solutionsequences.h: 
 */

#ifndef _SOLUTION_SEQUENCES_H_
#define _SOLUTION_SEQUENCES_H_

class FemModel;
class Parameters;
template <class doubletype> class Matrix;
template <class doubletype> class Vector;
#include "../shared/Numerics/types.h"

void solutionsequence_thermal_nonlinear(FemModel* femmodel);
void solutionsequence_hydro_nonlinear(FemModel* femmodel);
void solutionsequence_nonlinear(FemModel* femmodel,bool conserve_loads);
void solutionsequence_newton(FemModel* femmodel);
void solutionsequence_fct(FemModel* femmodel);
void solutionsequence_FScoupling_nonlinear(FemModel* femmodel,bool conserve_loads);
void solutionsequence_linear(FemModel* femmodel);
void solutionsequence_la(FemModel* femmodel);
void solutionsequence_la_theta(FemModel* femmodel);
void solutionsequence_adjoint_linear(FemModel* femmodel);

/*convergence*/
void convergence(bool* pconverged, Matrix<IssmDouble>* K_ff,Vector<IssmDouble>* p_f,Vector<IssmDouble>* u_f,Vector<IssmDouble>* u_f_old,IssmDouble eps_res,IssmDouble eps_rel,IssmDouble eps_abs);

#endif
