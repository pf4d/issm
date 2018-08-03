/*!\file: controlad_core.cpp
 * \brief: core of the ad control solution 
 */ 

#include <config.h>
#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

#if defined (_HAVE_M1QN3_)  & defined (_HAVE_ADOLC_)
/*m1qn3 prototypes {{{*/
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
void simulad(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs);
FemModel* presimulad(int* pintn, double** pX, FemModel* femmodel);
void postsimulad(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs);
/*}}}*/
void controlad_core(FemModel* femmodel){ /*{{{*/

	/*Intermediaries*/
	FemModel*    femmodelad=NULL;
	int          i;
	long         omode;
	IssmPDouble  f,dxmin,gttol;
	IssmDouble   dxmind,gttold; 
	int          maxsteps,maxiter;
	int          intn,solution_type;
	IssmPDouble  *X  = NULL;
	IssmDouble   *Xd  = NULL;
	IssmPDouble  *G  = NULL;

	/*Recover some parameters*/
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&maxsteps,InversionMaxstepsEnum);
	femmodel->parameters->FindParam(&maxiter,InversionMaxiterEnum);
	femmodel->parameters->FindParam(&dxmind,InversionDxminEnum); dxmin=reCast<IssmPDouble>(dxmind);
	femmodel->parameters->FindParam(&gttold,InversionGttolEnum); gttol=reCast<IssmPDouble>(gttold);
	femmodel->parameters->SetParam(false,SaveResultsEnum);

	/*Initialize M1QN3 parameters*/
	if(VerboseControl())_printf0_("   Initialize M1QN3 parameters\n");
	SimulFunc costfuncion  = &simulad;  /*Cost function address*/
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

	/*Run the first part of simulad, in order to get things started!:*/
	femmodelad=presimulad(&intn,&X,femmodel);

	/*Get problem dimension and initialize gradient: */
	long n = long(intn);
	G = xNew<IssmPDouble>(n);

	/*Allocate m1qn3 working arrays (see doc)*/
	long      m   = 100;
	long      ndz = 4*n+m*(2*n+1);
	double*   dz  = xNew<double>(ndz);

	if(VerboseControl())_printf0_("   Computing initial solution\n");
	_printf0_("\n");
	_printf0_("Cost function f(x)   | Gradient norm |g(x)| |  List of contributions\n");
	_printf0_("____________________________________________________________________\n");

	//run post simular phase, to fire up the control optimization
	postsimulad(&indic,&n,X,&f,G,izs,rzs,(void*)femmodelad); 
	double f1=f;

	m1qn3_(costfuncion,prosca,&ctonbe_,&ctcabe_,
				&n,X,&f,G,&dxmin,&f1,
				&gttol,normtype,&impres,&io,imode,&omode,&niter,&nsim,iz,dz,&ndz,
				&reverse,&indic,izs,rzs,(void*)femmodel);

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
	
	/*Save results:*/
	femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,AutodiffJacobianEnum,G,n,1,0,0.0));
	femmodel->results->AddObject(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,AutodiffXpEnum,X,intn,1,0,0.0));

	/*Clean-up and return*/
	xDelete<double>(G);
	xDelete<double>(X);
	xDelete<double>(dz);
}/*}}}*/
FemModel* presimulad(int* pintn, double** pX, FemModel* femmodel){ /*{{{*/

	/*Intermediaries:*/
	char* rootpath=NULL;
	char* inputfilename=NULL;
	char* outputfilename=NULL;
	char* toolkitsfilename=NULL;
	char* lockfilename=NULL;
	char* restartfilename=NULL;
	int         solution_type;
	IssmDouble    pfd;
	IssmDouble*   Xd=NULL;
	int           intn;
	IssmPDouble*   X=NULL;
	int            i;
	
	/*Now things get complicated. The femmodel we recovered did not initialize an AD trace, so we can't compute gradients with it. We are going to recreate 
	 *a new femmodel, identical in all aspects to the first one, with trace on though, which will allow us to run the forward mode and get the gradient 
	 in one run of the solution core. So first recover the filenames required for the FemModel constructor, then call a new ad tailored constructor:*/
	femmodel->parameters->FindParam(&rootpath,RootPathEnum);
	femmodel->parameters->FindParam(&inputfilename,InputFileNameEnum);
	femmodel->parameters->FindParam(&outputfilename,OutputFileNameEnum);
	femmodel->parameters->FindParam(&toolkitsfilename,ToolkitsFileNameEnum);
	femmodel->parameters->FindParam(&lockfilename,LockFileNameEnum);
	femmodel->parameters->FindParam(&restartfilename,RestartFileNameEnum);

	femmodel=new FemModel(rootpath, inputfilename, outputfilename, toolkitsfilename, lockfilename, restartfilename,IssmComm::GetComm(), femmodel->solution_type,NULL);

	
	/*Get initial guess:*/
	femmodel->parameters->FindParam(&Xd,&intn,AutodiffXpEnum);
	X=xNew<IssmPDouble>(intn); for(i=0;i<intn;i++) X[i]=reCast<IssmPDouble>(Xd[i]); 

	xDelete<char>(rootpath);
	xDelete<char>(inputfilename);
	xDelete<char>(outputfilename);
	xDelete<char>(toolkitsfilename);
	xDelete<char>(lockfilename);
	xDelete<char>(restartfilename);
	xDelete<IssmDouble>(Xd);

	*pintn=intn;
	*pX=X;

	return femmodel;

} /*}}}*/
void postsimulad(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs){ /*{{{*/

	/*Intermediaries:*/
	char* rootpath=NULL;
	char* inputfilename=NULL;
	char* outputfilename=NULL;
	char* toolkitsfilename=NULL;
	char* lockfilename=NULL;
	IssmPDouble* G2=NULL;
	int         solution_type;
	FemModel   *femmodel  =  NULL;
	IssmDouble    pfd;
	int            i;
	
	/*Recover Femmodel*/
	femmodel  = (FemModel*)dzs;

	/*Recover number of cost functions responses*/
	int num_responses;
	femmodel->parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);

	/*Recover some parameters*/
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);

	/*Compute solution:*/
	void (*solutioncore)(FemModel*)=NULL;
	CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);
	solutioncore(femmodel);

	/*Compute objective function*/
	IssmDouble* Jlist = NULL;
	femmodel->CostFunctionx(&pfd,&Jlist,NULL); *pf=reCast<IssmPDouble>(pfd);
	_printf0_("f(x) = "<<setw(12)<<setprecision(7)<<*pf<<"  |  ");
	
	/*Compute gradient using AD. Gradient is in the results after the ad_core is called*/
	adgradient_core(femmodel); 

	if(IssmComm::GetRank()==0){
		IssmPDouble* G_temp=NULL;
		GenericExternalResult<IssmPDouble*>* gradient=(GenericExternalResult<IssmPDouble*>*)femmodel->results->FindResult(AutodiffJacobianEnum); _assert_(gradient);
		G_temp=gradient->GetValues();
		/*copy onto G2, to avoid a leak: */
		G2=xNew<IssmPDouble>(*n); 
		xMemCpy<IssmPDouble>(G2,G_temp,*n);
	}
	else G2=xNew<IssmPDouble>(*n);

	/*MPI broadcast results:*/
	ISSM_MPI_Bcast(G2,*n,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
	
	/*Send gradient to m1qn3 core: */
	for(long i=0;i<*n;i++) G[i] = G2[i];
	
	/*Constrain X and G*/
	IssmDouble  Gnorm = 0.;
	for(long i=0;i<*n;i++) Gnorm += G[i]*G[i];
	Gnorm = sqrt(Gnorm);

	/*Print info*/
	_printf0_("       "<<setw(12)<<setprecision(7)<<Gnorm<<" |");
	for(int i=0;i<num_responses;i++) _printf0_(" "<<setw(12)<<setprecision(7)<<Jlist[i]);
	_printf0_("\n");

	/*Clean-up and return*/
	xDelete<IssmDouble>(Jlist);
	xDelete<IssmPDouble>(G2);
	
	xDelete<char>(rootpath);
	xDelete<char>(inputfilename);
	xDelete<char>(outputfilename);
	xDelete<char>(toolkitsfilename);
	xDelete<char>(lockfilename);

} /*}}}*/
void simulad(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs){ /*{{{*/

	/*Intermediaries:*/
	char* rootpath=NULL;
	char* inputfilename=NULL;
	char* outputfilename=NULL;
	char* toolkitsfilename=NULL;
	char* lockfilename=NULL;
	char* restartfilename=NULL;
	IssmPDouble* G2=NULL;
	int         solution_type;
	FemModel   *femmodel  =  NULL;
	FemModel   *femmodelad  = NULL;
	IssmDouble    pfd;
	int            i;
	
	/*Recover Femmodel*/
	femmodel  = (FemModel*)dzs;

	/*Recover number of cost functions responses*/
	int num_responses;
	femmodel->parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);

	/*Now things get complicated. The femmodel we recovered did not initialize an AD trace, so we can't compute gradients with it. We are going to recreate 
	 *a new femmodel, identical in all aspects to the first one, with trace on though, which will allow us to run the forward mode and get the gradient 
	 in one run of the solution core. So first recover the filenames required for the FemModel constructor, then call a new ad tailored constructor:*/
	femmodel->parameters->FindParam(&rootpath,RootPathEnum);
	femmodel->parameters->FindParam(&inputfilename,InputFileNameEnum);
	femmodel->parameters->FindParam(&outputfilename,OutputFileNameEnum);
	femmodel->parameters->FindParam(&toolkitsfilename,ToolkitsFileNameEnum);
	femmodel->parameters->FindParam(&lockfilename,LockFileNameEnum);
	femmodel->parameters->FindParam(&restartfilename,RestartFileNameEnum);

	femmodelad=new FemModel(rootpath, inputfilename, outputfilename, toolkitsfilename, lockfilename, restartfilename,IssmComm::GetComm(), femmodel->solution_type,X);
	femmodel=femmodelad; //We can do this, because femmodel is being called from outside, not by reference, so we won't erase it
	
	/*Recover some parameters*/
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);

	/*Compute solution:*/
	void (*solutioncore)(FemModel*)=NULL;
	CorePointerFromSolutionEnum(&solutioncore,femmodel->parameters,solution_type);
	solutioncore(femmodel);

	/*Compute objective function*/
	IssmDouble* Jlist = NULL;
	femmodel->CostFunctionx(&pfd,&Jlist,NULL); *pf=reCast<IssmPDouble>(pfd);
	_printf0_("f(x) = "<<setw(12)<<setprecision(7)<<*pf<<"  |  ");
	
	/*Compute gradient using AD. Gradient is in the results after the ad_core is called*/
	adgradient_core(femmodel); 

	if(IssmComm::GetRank()==0){
		IssmPDouble* G_temp=NULL;
		GenericExternalResult<IssmPDouble*>* gradient=(GenericExternalResult<IssmPDouble*>*)femmodel->results->FindResult(AutodiffJacobianEnum); _assert_(gradient);
		G_temp=gradient->GetValues();
		/*copy onto G2, to avoid a leak: */
		G2=xNew<IssmPDouble>(*n); 
		xMemCpy<IssmPDouble>(G2,G_temp,*n);
	}
	else G2=xNew<IssmPDouble>(*n);

	/*MPI broadcast results:*/
	ISSM_MPI_Bcast(G2,*n,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
	
	/*Send gradient to m1qn3 core: */
	for(long i=0;i<*n;i++) G[i] = G2[i];
	
	/*Recover Gnorm: */
	IssmDouble  Gnorm = 0.;
	for(int i=0;i<*n;i++) Gnorm += G[i]*G[i];
	Gnorm = sqrt(Gnorm);

	/*Print info*/
	_printf0_("       "<<setw(12)<<setprecision(7)<<Gnorm<<" |");
	for(int i=0;i<num_responses;i++) _printf0_(" "<<setw(12)<<setprecision(7)<<Jlist[i]);
	_printf0_("\n");

	/*Clean-up and return*/
	xDelete<IssmDouble>(Jlist);
	xDelete<IssmPDouble>(G2);
	
	xDelete<char>(rootpath);
	xDelete<char>(inputfilename);
	xDelete<char>(outputfilename);
	xDelete<char>(toolkitsfilename);
	xDelete<char>(lockfilename);
	xDelete<char>(restartfilename);
	if(femmodelad)delete femmodelad;

} /*}}}*/
#else
void controlad_core(FemModel* femmodel){ /*{{{*/
	_error_("AD and/or M1QN3 not installed");
}/*}}}*/
#endif //_HAVE_M1QN3_
