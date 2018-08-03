/*!\file: sealevelrise_core_noneustatic.cpp //this computes the contributions from Eq.4 of Farrel and Clarke, rhs terms 2 and 5.
 * \brief: non eustatic core of the SLR solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void slrconvergence(bool* pconverged, Vector<IssmDouble>* Sg,Vector<IssmDouble>* Sg_old,IssmDouble eps_rel,IssmDouble eps_abs);

Vector<IssmDouble>* sealevelrise_core_noneustatic(FemModel* femmodel,Vector<IssmDouble>* Sg_eustatic){ /*{{{*/

	Vector<IssmDouble> *Sg    = NULL;
	Vector<IssmDouble> *Sg_old    = NULL;

	Vector<IssmDouble> *Sgo    = NULL; //ocean convolution of the perturbation to gravity potential.
	Vector<IssmDouble> *Sgo_rot= NULL; // rotational feedback 
	IssmDouble          Sgo_oceanaverage = 0;  //average of Sgo over the ocean.

	/*parameters: */
	int count;
	bool save_results;
	int  gsize;
	int  configuration_type;
	bool spherical=true;
	bool converged=true;
	bool rotation=true;
	bool verboseconvolution=true;
	int max_nonlinear_iterations;
	IssmDouble           eps_rel;
	IssmDouble           eps_abs;
	IssmDouble          *latitude    = NULL;
	IssmDouble          *longitude    = NULL;
	IssmDouble          *radius    = NULL;
	IssmDouble           eustatic;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&max_nonlinear_iterations,SealevelriseMaxiterEnum);
	femmodel->parameters->FindParam(&eps_rel,SealevelriseReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,SealevelriseAbstolEnum);
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	
	/*computational flag: */
	femmodel->parameters->FindParam(&rotation,SealevelriseRotationEnum);

	/*first, recover lat,long and radius vectors from vertices: */
	VertexCoordinatesx(&latitude,&longitude,&radius,femmodel->vertices,spherical); 

	/*Figure out size of g-set deflection vector and allocate solution vector: */
	gsize = femmodel->nodes->NumberOfDofs(configuration_type,GsetEnum);
	
	/*Initialize:*/
	Sg = new Vector<IssmDouble>(gsize);
	Sg->Assemble();
	Sg_eustatic->Copy(Sg);  //first initialize Sg with the eustatic component computed in sealevelrise_core_eustatic.

	Sg_old = new Vector<IssmDouble>(gsize);
	Sg_old->Assemble();

	count=1;
	converged=false;
	
	/*Start loop: */
	for(;;){

		//save pointer to old sea level rise
		delete Sg_old; Sg_old=Sg; 

		/*Initialize solution vector: */
		Sg  = new Vector<IssmDouble>(gsize); Sg->Assemble();
		Sgo = new Vector<IssmDouble>(gsize); Sgo->Assemble();

		/*call the non eustatic module: */
		femmodel->SealevelriseNonEustatic(Sgo, Sg_old, latitude, longitude, radius,verboseconvolution);
	
		/*assemble solution vector: */
		Sgo->Assemble(); 

		if(rotation){
			/*call rotational feedback  module: */
			Sgo_rot = new Vector<IssmDouble>(gsize); Sgo_rot->Assemble();
			femmodel->SealevelriseRotationalFeedback(Sgo_rot,Sg_old,latitude,longitude,radius); 
			Sgo_rot->Assemble(); 

			Sgo->AXPY(Sgo_rot,1); 
		}
		
		/*we need to average Sgo over the ocean: RHS term  5 in Eq.4 of Farrel and clarke. Only the elements can do that: */
		Sgo_oceanaverage=femmodel->SealevelriseOceanAverage(Sgo);
	
		/*Sg is the sum of the eustatic term, and the ocean terms: */
		Sg_eustatic->Copy(Sg); Sg->AXPY(Sgo,1); 
		Sg->Shift(-Sgo_oceanaverage);

		/*convergence criterion:*/
		slrconvergence(&converged,Sg,Sg_old,eps_rel,eps_abs);

		/*Increase count: */
		count++;
		if(converged==true){
			break;
		}
		if(count>=max_nonlinear_iterations){
			_printf0_("   maximum number of nonlinear iterations (" << max_nonlinear_iterations << ") exceeded\n"); 
			converged=true;
			break;
		}	
		
		/*some minor verbosing adjustment:*/
		if(count>1)verboseconvolution=false;
		
	}
	if(VerboseConvergence()) _printf0_("\n   total number of iterations: " << count-1 << "\n");

	xDelete<IssmDouble>(latitude);
	xDelete<IssmDouble>(longitude);
	xDelete<IssmDouble>(radius);
	delete Sg_old;

	return Sg;
} /*}}}*/

void slrconvergence(bool* pconverged, Vector<IssmDouble>* Sg,Vector<IssmDouble>* Sg_old,IssmDouble eps_rel,IssmDouble eps_abs){ /*{{{*/
	
	bool converged=true;
	IssmDouble ndS,nS; 
	Vector<IssmDouble> *dSg    = NULL;

	//compute norm(du) and norm(u) if requested
	dSg=Sg_old->Duplicate(); Sg_old->Copy(dSg); dSg->AYPX(Sg,-1.0);
	ndS=dSg->Norm(NORM_TWO); 
	
	if(!xIsNan<IssmDouble>(eps_rel)){
		nS=Sg_old->Norm(NORM_TWO);
	}

	if (xIsNan<IssmDouble>(ndS) || xIsNan<IssmDouble>(nS)) _error_("convergence criterion is NaN!");

	//clean up
	delete dSg;

	//print
	if(!xIsNan<IssmDouble>(eps_rel)){
		if((ndS/nS)<eps_rel){
			if(VerboseConvergence()) _printf0_(setw(50) << left << "      convergence criterion: norm(dS)/norm(S)" << ndS/nS*100 << " < " << eps_rel*100 << " %\n");
		}
		else{ 
			if(VerboseConvergence()) _printf0_(setw(50) << left << "      convergence criterion: norm(dS)/norm(S)" << ndS/nS*100 << " > " << eps_rel*100 << " %\n");
			converged=false;
		}
	}
	if(!xIsNan<IssmDouble>(eps_abs)){
		if(ndS<eps_abs){
			if(VerboseConvergence()) _printf0_(setw(50) << left << "      convergence criterion: norm(dS)" << ndS << " < " << eps_abs << " \n");
		}
		else{ 
			if(VerboseConvergence()) _printf0_(setw(50) << left << "      convergence criterion: norm(dS)" << ndS << " > " << eps_abs << " \n");
			converged=false;
		}
	}

	/*assign output*/
	*pconverged=converged;

} /*}}}*/

