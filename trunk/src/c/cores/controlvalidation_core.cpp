/*!\file: controlvalidation_core.cpp
 * \brief: core of the control solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void controlvalidation_core(FemModel* femmodel){

	int         solution_type,n;
	int         num_responses;
	IssmDouble  j0,j,yts;
	IssmDouble  Ialpha,exponent,alpha;
	IssmDouble* scaling_factors = NULL;
	IssmDouble* jlist = NULL;
	IssmDouble *G = NULL;
	IssmDouble *X = NULL;
	IssmDouble *X0= NULL;

	/*Solution and Adjoint core pointer*/
	void (*solutioncore)(FemModel*) = NULL;
	void (*adjointcore)(FemModel*)  = NULL;

	/*Recover parameters used throughout the solution*/
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->SetParam(false,SaveResultsEnum);
	femmodel->parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);
	femmodel->parameters->FindParam(&scaling_factors,NULL,InversionControlScalingFactorsEnum);

	/*Get initial guess*/
	Vector<IssmDouble> *Xpetsc = NULL;
	GetVectorFromControlInputsx(&Xpetsc,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"value");
	Xpetsc->GetSize(&n);
	X0 = Xpetsc->ToMPISerial();
	delete Xpetsc;

	/*Allocate current vector*/
	X = xNew<IssmDouble>(n);

	/*out of solution_type, figure out solution core and adjoint function pointer*/
	CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);
	AdjointCorePointerFromSolutionEnum(&adjointcore,solution_type);

	if(VerboseControl()) _printf0_("   Compute Initial solution\n");
	solutioncore(femmodel);
	if(VerboseControl()) _printf0_("   Compute Adjoint\n");
	adjointcore(femmodel);

	if(VerboseControl()) _printf0_("   Compute Initial cost function\n");
	femmodel->CostFunctionx(&j0,&jlist,NULL);
	_printf0_("Initial cost function J(x) = "<<setw(12)<<setprecision(7)<<j0<<"\n");
	xDelete<IssmDouble>(jlist);

	if(VerboseControl()) _printf0_("   Compute Gradient\n");
	Gradjx(&G,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters);
	for(int i=0;i<n;i++) G[i] = -G[i];

	/*Allocate output*/
	int num = 26;
	IssmDouble* output = xNew<IssmDouble>(2*num);

	/*Start loop*/
	_printf0_("       alpha      Ialpha \n");
	_printf0_("_________________________\n");
	for(int m=0;m<num;m++){

		/*Create new vector*/
		alpha    = pow(2.,-m);
		for(int i=0;i<n;i++) X[i] = X0[i] + alpha*scaling_factors[0];

		/*Calculate j(k+alpha delta k) */
		SetControlInputsFromVectorx(femmodel,X);
		solutioncore(femmodel);
		femmodel->CostFunctionx(&j,NULL,NULL);

		IssmDouble Den = 0.;
		for(int i=0;i<n;i++) Den += alpha* G[i] * scaling_factors[0];
		Ialpha = fabs((j - j0)/Den - 1.);

		_printf0_(" " << setw(11) << setprecision (5)<<alpha<<" " << setw(11) << setprecision (5)<<Ialpha<<"\n");
		output[m*2+0] = alpha;
		output[m*2+1] = Ialpha;
	}

	/*output*/
	#ifdef _HAVE_ADOLC_
	IssmPDouble* J_passive=xNew<IssmPDouble>(2*num);
	for(int i=0;i<2*num;i++) J_passive[i]=reCast<IssmPDouble>(output[i]);
	femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,JEnum,J_passive,num,2,0,0));
	xDelete<IssmPDouble>(J_passive);
	#else
	femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,JEnum,output,num,2,0,0));
	#endif
	ControlInputSetGradientx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,G);
	femmodel->OutputControlsx(&femmodel->results);

	/*Clean up and return*/
	xDelete<IssmDouble>(output);
	xDelete<IssmDouble>(G);
	xDelete<IssmDouble>(X);
	xDelete<IssmDouble>(X0);
}
