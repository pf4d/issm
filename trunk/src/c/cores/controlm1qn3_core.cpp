/*!\file: controlm1qn3_core.cpp
 * \brief: core of the control solution 
 */ 

#include <config.h>
#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

#if defined (_HAVE_M1QN3_) & !defined(_HAVE_ADOLC_)
/*m1qn3 prototypes*/
extern "C" void *ctonbe_; // DIS mode : Conversion
extern "C" void *ctcabe_; // DIS mode : Conversion
extern "C" void *euclid_; // Scalar product
typedef void (*SimulFunc) (long* indic,long* n, double* x, double* pf,double* g,long [],float [],void* dzs);
extern "C" void m1qn3_ (void f(long* indic,long* n, double* x, double* pf,double* g,long [],float [],void* dzs),
			void **, void **, void **,
			long *, double [], double *, double [], double*, double *,
			double *, char [], long *, long *, long *, long *, long *, long *, long [], double [], long *,
			long *, long *, long [], float [],void* );

/*Cost function prototype*/
void simul(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs);

/*Use struct to provide arguments*/
typedef struct{
	FemModel  * femmodel;
	IssmDouble* Jlist;
	int         M;
	int         N;
	int*        i;
} m1qn3_struct;

void controlm1qn3_core(FemModel* femmodel){

	/*Intermediaries*/
	long         omode;
	double       f,dxmin,gttol; 
	int          maxsteps,maxiter;
	int          intn,numberofvertices,num_controls,num_cost_functions,solution_type;
	IssmDouble  *scaling_factors = NULL;
	IssmDouble  *X  = NULL;
	IssmDouble  *G  = NULL;

	/*Recover some parameters*/
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	femmodel->parameters->FindParam(&num_cost_functions,InversionNumCostFunctionsEnum);
	femmodel->parameters->FindParam(&maxsteps,InversionMaxstepsEnum);
	femmodel->parameters->FindParam(&maxiter,InversionMaxiterEnum);
	femmodel->parameters->FindParam(&dxmin,InversionDxminEnum);
	femmodel->parameters->FindParam(&gttol,InversionGttolEnum);
	femmodel->parameters->FindParam(&scaling_factors,NULL,InversionControlScalingFactorsEnum);
	femmodel->parameters->SetParam(false,SaveResultsEnum);
	numberofvertices=femmodel->vertices->NumberOfVertices();

	/*Initialize M1QN3 parameters*/
	if(VerboseControl())_printf0_("   Initialize M1QN3 parameters\n");
	SimulFunc costfuncion  = &simul;    /*Cost function address*/
	void**    prosca       = &euclid_;  /*Dot product function (euclid is the default)*/
	char      normtype[]   = "dfn";     /*Norm type: dfn = scalar product defined by prosca*/
	long      izs[5];                   /*Arrays used by m1qn3 subroutines*/
	long      iz[5];                    /*Integer m1qn3 working array of size 5*/
	float     rzs[1];                   /*Arrays used by m1qn3 subroutines*/
	long      impres       = 0;         /*verbosity level*/
	long      imode[3]     = {0};       /*scaling and starting mode, 0 by default*/
	long      indic        = 4;         /*compute f and g*/
	long      reverse      = 0;         /*reverse or direct mode*/
	long      io           = 6;         /*Channel number for the output*/

	/*Optimization criterions*/
	long niter = long(maxsteps); /*Maximum number of iterations*/
	long nsim  = long(maxiter);/*Maximum number of function calls*/

	/*Get initial guess*/
	Vector<IssmDouble> *Xpetsc = NULL;
	GetVectorFromControlInputsx(&Xpetsc,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"value");
	X = Xpetsc->ToMPISerial();
	Xpetsc->GetSize(&intn);
	delete Xpetsc;
	_assert_(intn==numberofvertices*num_controls);

	/*Get problem dimension and initialize gradient and initial guess*/
	long n = long(intn);
	G = xNew<double>(n);

	/*Scale control for M1QN3*/
	for(int i=0;i<numberofvertices;i++){
		for(int c=0;c<num_controls;c++){
			int index = num_controls*i+c;
			X[index] = X[index]/scaling_factors[c];
		}
	}

	/*Allocate m1qn3 working arrays (see documentation)*/
	long      m   = 100;
	long      ndz = 4*n+m*(2*n+1);
	double*   dz  = xNew<double>(ndz);

	if(VerboseControl())_printf0_("   Computing initial solution\n");
	_printf0_("\n");
	_printf0_("Cost function f(x)   | Gradient norm |g(x)| |  List of contributions\n");
	_printf0_("____________________________________________________________________\n");

	/*Prepare structure for m1qn3*/
	m1qn3_struct mystruct;
	mystruct.femmodel = femmodel;
	mystruct.M        = maxiter;
	mystruct.N        = num_cost_functions+1;
	mystruct.Jlist    = xNewZeroInit<IssmDouble>(mystruct.M*mystruct.N);
	mystruct.i        = xNewZeroInit<int>(1);

	/*Initialize Gradient and cost function of M1QN3*/
	indic = 4; /*gradient required*/
	simul(&indic,&n,X,&f,G,izs,rzs,(void*)&mystruct);

	/*Estimation of the expected decrease in f during the first iteration*/
	double df1=f;

	/*Call M1QN3 solver*/
	m1qn3_(costfuncion,prosca,&ctonbe_,&ctcabe_,
				&n,X,&f,G,&dxmin,&df1,
				&gttol,normtype,&impres,&io,imode,&omode,&niter,&nsim,iz,dz,&ndz,
				&reverse,&indic,izs,rzs,(void*)&mystruct);

	switch(int(omode)){
		case 0:  _printf0_("   Stop requested (indic = 0)\n"); break;
		case 1:  _printf0_("   Convergence reached (gradient satisfies stopping criterion)\n"); break;
		case 2:  _printf0_("   Bad initialization\n"); break;
		case 3:  _printf0_("   Line search failure\n"); break;
		case 4:  _printf0_("   Maximum number of iterations exceeded\n");break;
		case 5:  _printf0_("   Maximum number of function calls exceeded\n"); break;
		case 6:  _printf0_("   stopped on dxmin during line search\n"); break;
		case 7:  _printf0_("   <g,d> > 0  or  <y,s> <0\n"); break;
		default: _printf0_("   Unknown end condition\n");
	}

	/*Constrain solution vector*/
	IssmDouble  *XL = NULL;
	IssmDouble  *XU = NULL;
	GetVectorFromControlInputsx(&XL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"lowerbound");
	GetVectorFromControlInputsx(&XU,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"upperbound");
	for(int i=0;i<numberofvertices;i++){
		for(int c=0;c<num_controls;c++){
			int index = num_controls*i+c;
			X[index] = X[index]*scaling_factors[c];
			if(X[index]>XU[index]) X[index]=XU[index];
			if(X[index]<XL[index]) X[index]=XL[index];
		}
	}
	SetControlInputsFromVectorx(femmodel,X);
	ControlInputSetGradientx(femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,G);
	femmodel->OutputControlsx(&femmodel->results);
	femmodel->results->AddObject(new GenericExternalResult<double*>(femmodel->results->Size()+1,JEnum,mystruct.Jlist,(*mystruct.i),mystruct.N,0,0));

	/*Finalize*/
	if(VerboseControl()) _printf0_("   preparing final solution\n");
	femmodel->parameters->SetParam(true,SaveResultsEnum);
	void (*solutioncore)(FemModel*)=NULL;
	CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);
	solutioncore(femmodel);

	/*Clean-up and return*/
	xDelete<double>(G);
	xDelete<double>(X);
	xDelete<double>(dz);
	xDelete<double>(XU);
	xDelete<double>(XL);
	xDelete<double>(scaling_factors);
	xDelete<double>(mystruct.Jlist);
	xDelete<int>(mystruct.i);
}

/*Cost function definition*/
void simul(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs){

	/*Recover Arguments*/
	m1qn3_struct *input_struct = (m1qn3_struct*)dzs;
	FemModel     *femmodel     = input_struct->femmodel;
	IssmDouble   *Jlist        = input_struct->Jlist;
	int           JlistM       = input_struct->M;
	int           JlistN       = input_struct->N;
	int          *Jlisti       = input_struct->i;

	/*Recover some parameters*/
	int num_responses,num_controls,numberofvertices,solution_type;
	IssmDouble* scaling_factors = NULL;
	femmodel->parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	femmodel->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	femmodel->parameters->FindParam(&scaling_factors,NULL,InversionControlScalingFactorsEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	numberofvertices=femmodel->vertices->NumberOfVertices();

	/*Constrain input vector and update controls*/
	IssmDouble  *XL = NULL;
	IssmDouble  *XU = NULL;
	GetVectorFromControlInputsx(&XL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"lowerbound");
	GetVectorFromControlInputsx(&XU,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"upperbound");
	for(int i=0;i<numberofvertices;i++){
		for(int c=0;c<num_controls;c++){
			int index = num_controls*i+c;
			X[index] = X[index]*scaling_factors[c];
			if(X[index]>XU[index]) X[index]=XU[index];
			if(X[index]<XL[index]) X[index]=XL[index];
		}
	}
	SetControlInputsFromVectorx(femmodel,X);

	/*Compute solution and adjoint*/
	void (*solutioncore)(FemModel*)=NULL;
	void (*adjointcore)(FemModel*)=NULL;
	CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);
	solutioncore(femmodel);

	/*Check size of Jlist to avoid crashes*/
	_assert_((*Jlisti)<JlistM);
	_assert_(JlistN==num_responses+1);

	/*Compute objective function*/
	IssmDouble* Jtemp = NULL;
	femmodel->CostFunctionx(pf,&Jtemp,NULL);
	_printf0_("f(x) = "<<setw(12)<<setprecision(7)<<*pf<<"  |  ");

	/*Record cost function values and delete Jtemp*/
	for(int i=0;i<num_responses;i++) Jlist[(*Jlisti)*JlistN+i] = Jtemp[i];
	Jlist[(*Jlisti)*JlistN+num_responses] = *pf;
	xDelete<IssmDouble>(Jtemp);

	if(*indic==0){
		/*dry run, no gradient required*/

		/*Retrieve objective functions independently*/
		_printf0_("            N/A |\n");
		for(int i=0;i<num_responses;i++) _printf0_(" "<<setw(12)<<setprecision(7)<<Jlist[(*Jlisti)*JlistN+i]);
		_printf0_("\n");

		*Jlisti = (*Jlisti) +1;
		xDelete<IssmDouble>(XU);
		xDelete<IssmDouble>(XL);
		return;
	}

	/*Compute Adjoint*/
	AdjointCorePointerFromSolutionEnum(&adjointcore,solution_type);
	adjointcore(femmodel);

	/*Compute gradient*/
	IssmDouble* G2 = NULL;
	Gradjx(&G2,NULL,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters);
	for(long i=0;i<*n;i++) G[i] = -G2[i];
	xDelete<IssmDouble>(G2);

	/*Constrain Gradient*/
	IssmDouble  Gnorm = 0.;
	for(int i=0;i<numberofvertices;i++){
		for(int c=0;c<num_controls;c++){
			int index = num_controls*i+c;
			if(X[index]>=XU[index]) G[index]=0.;
			if(X[index]<=XL[index]) G[index]=0.;
			G[index] = G[index]*scaling_factors[c];
			X[index] = X[index]/scaling_factors[c];
			Gnorm += G[index]*G[index];
		}
	}
	Gnorm = sqrt(Gnorm);

	/*Print info*/
	_printf0_("       "<<setw(12)<<setprecision(7)<<Gnorm<<" |");
	for(int i=0;i<num_responses;i++) _printf0_(" "<<setw(12)<<setprecision(7)<<Jlist[(*Jlisti)*JlistN+i]);
	_printf0_("\n");

	/*Clean-up and return*/
	*Jlisti = (*Jlisti) +1;
	xDelete<IssmDouble>(XU);
	xDelete<IssmDouble>(XL);
	xDelete<IssmDouble>(scaling_factors);
}

#else
void controlm1qn3_core(FemModel* femmodel){
	_error_("M1QN3 not installed");
}
#endif //_HAVE_M1QN3_
