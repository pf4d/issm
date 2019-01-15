/*!\file: steadystate_core.cpp
 * \brief: core of the steadystate solution
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

// FIXME:  use a header.
// local prototype :
IssmDouble calc_relative_error(Vector<IssmDouble>* u,
                               Vector<IssmDouble>* u_old);

bool steadystateconvergence(Vector<IssmDouble>* tg,
                            Vector<IssmDouble>* tg_old,
                            Vector<IssmDouble>* ug,
                            Vector<IssmDouble>* ug_old,
                            IssmDouble          reltol);

void steadystate_core(FemModel* femmodel)
{
	/*intermediary: */
	Vector<IssmDouble>* ug     = NULL;
	Vector<IssmDouble>* ug_old = NULL;
	Vector<IssmDouble>* tg     = NULL;
	Vector<IssmDouble>* tg_old = NULL;
	Vector<IssmDouble>* tt     = NULL;
	Vector<IssmDouble>* tt_old = NULL;

	/*parameters: */
	bool        save_results, isenthalpy;
	bool        converged, thermal_converged;
	int         maxiter;
	int         thermal_step;
	IssmDouble  reltol;
	IssmDouble  t_rel_err;
	IssmDouble  t_reltol = 1e-15;
	int         step;
	int         numoutputs        = 0;
	char**      requested_outputs = NULL;

	/* recover parameters:*/
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&maxiter,SteadystateMaxiterEnum);
	femmodel->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);
	femmodel->parameters->FindParam(&reltol,SteadystateReltolEnum);
	femmodel->parameters->SetParam(false,SaveResultsEnum);
	femmodel->parameters->FindParam(&numoutputs,SteadystateNumRequestedOutputsEnum);

	if(numoutputs)
	{
		femmodel->parameters->FindParam(&requested_outputs,&numoutputs,SteadystateRequestedOutputsEnum);
	}

	// begin the algorithm :
	converged = false;
	step      = 1;
	while(!converged && step < maxiter)
	{
		// Compute first velocity, then temperature due to high
		// sensitivity of temperature to velocity :
		if(VerboseSolution())
		{
			_printf0_("\n");
			_printf0_("======================================================\n");
			_printf0_("   computing velocity and temperature for step: " << step << "\n");
			_printf0_("======================================================\n");
		}

		if(VerboseSolution()) _printf0_("\n   -- computing new velocity --\n\n");
		stressbalance_core(femmodel);
		GetSolutionFromInputsx(&ug,femmodel);

		if(VerboseSolution()) _printf0_("\n   -- computing new temperature --\n\n");

		// quasi-time step until the thermal solution has converged :
		thermal_converged = false;
		thermal_step      = 1;
		while (!thermal_converged)
		{
			thermal_core(femmodel);
			if (!isenthalpy)
			{
				// could also be MeltingAnalysis :
				femmodel->SetCurrentConfiguration(ThermalAnalysisEnum);
			}
			GetSolutionFromInputsx(&tt,femmodel);

			if (thermal_step > 1)
			{
				// compute norm(dt) / norm(t) :
				t_rel_err = calc_relative_error(tt, tt_old);
				if (xIsNan<IssmDouble>(t_rel_err))
				{
					_error_("convergence criterion is NaN!");
				}
				else if ( t_rel_err < t_reltol )
				{
					if (VerboseConvergence())
					{
						_printf0_(setw(50)
						          << left
						          << "\n   Temperature convergence: ||dT|| / ||T||"
						          << t_rel_err*100
						          << " < "
						          << t_reltol*100
						          << " %\n\n");
					}
					thermal_converged = true;
				}
				else
				{
					if (VerboseConvergence())
					{
						_printf0_(setw(50)
						          << left
						          << "\n   Temperature convergence: ||dT|| / ||T||"
						          << t_rel_err*100
						          << " > "
						          << t_reltol*100
						          << " %\n\n");
					}
				}
			}
			//delete tt_old; tt_old=tt;
			tt_old = tt->Duplicate();
			tt->Copy(tt_old);
			thermal_step++;
		}
		tg = tt->Duplicate();
		tt->Copy(tg);
		if (step > 1)
		{
			converged = steadystateconvergence(tg, tg_old, ug, ug_old, reltol);
		}
		if (step > (maxiter - 1))
		{
			if(VerboseSolution())
			{
				_printf0_("   maximum number steadystate iterations " << maxiter << " reached\n");
			}
		}
		// update results and increase counter :
		//delete tg_old; tg_old=tg;
		//delete ug_old; ug_old=ug;
		tg_old = tg->Duplicate();
		tg->Copy(tg_old);
		ug_old = ug->Duplicate();
		ug->Copy(ug_old);
		step++;
	}

	if(save_results)
	{
		if(VerboseSolution())
		{
			_printf0_("   saving results\n");
		}
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
	}

	// free resources :
	delete tt_old;
	delete tg_old;
	delete ug_old;
	delete tt;
	delete tg;
	delete ug;

	if(numoutputs)
	{
		for(int i=0; i < numoutputs; i++)
		{
			xDelete<char>(requested_outputs[i]);
		}
		xDelete<char*>(requested_outputs);
	}
}

IssmDouble calc_relative_error(Vector<IssmDouble>* u,
                               Vector<IssmDouble>* u_old)
{
	// intermediary variables :
	Vector<IssmDouble>* du    = NULL;
	IssmDouble          ndu,nu;

	// compute norm(du) / norm(u) :
	du = u_old->Duplicate();
	u_old->Copy(du);
	du->AYPX(u,-1.0);
	ndu = du->Norm(NORM_TWO);
	nu  = u_old->Norm(NORM_TWO);
	return (ndu / nu);
}

bool steadystateconvergence(Vector<IssmDouble>* tg,
                            Vector<IssmDouble>* tg_old,
                            Vector<IssmDouble>* ug,
                            Vector<IssmDouble>* ug_old,
                            IssmDouble          reltol)
{
	// output :
	bool converged = true;

	// intermediary variables :
	Vector<IssmDouble>* dug    = NULL;
	Vector<IssmDouble>* dtg    = NULL;
	IssmDouble          u_rel_err, t_rel_err;

	if(VerboseSolution())
	{
		_printf0_("   checking steadystate convergence\n");
	}

	// compute norm(du) / norm(u) :
	u_rel_err = calc_relative_error(ug, ug_old);
	if (xIsNan<IssmDouble>(u_rel_err))
	{
		_error_("convergence criterion is NaN!");
	}
	else if ( u_rel_err < reltol )
	{
		if (VerboseConvergence())
		{
			_printf0_("\n"
			          << setw(50)
			          << left
			          << "   Velocity convergence:    ||du|| / ||u||"
			          << u_rel_err*100
			          << " < "
			          << reltol*100
			          << " %\n");
		}
	}
	else
	{
		if (VerboseConvergence())
		{
			_printf0_("\n"
			          << setw(50)
			          << left
			          << "   Velocity convergence:    ||du|| / ||u||"
			          << u_rel_err*100
			          << " > "
			          << reltol*100
			          << " %\n");
		}
		converged = false;
	}

	// compute norm(dt) / norm(t) :
	t_rel_err = calc_relative_error(tg, tg_old);
	if (xIsNan<IssmDouble>(t_rel_err))
	{
		_error_("convergence criterion is NaN!");
	}
	else if ( t_rel_err < reltol )
	{
		if(VerboseConvergence())
		{
			_printf0_(setw(50)
			          << left
			          << "   Temperature convergence: ||dT|| / ||T||"
			          << t_rel_err*100
			          << " < "
			          << reltol*100
			          << " %\n");
		}
	}
	else
	{
		if (VerboseConvergence())
		{
			_printf0_(setw(50)
			          << left
			          << "   Temperature convergence: ||dT|| / ||T||"
			          << t_rel_err*100
			          << " > "
			          << reltol*100
			          << " %\n");
		}
		converged = false;
	}

	// clean up and return
	delete dtg;
	delete dug;
	return converged;
}



