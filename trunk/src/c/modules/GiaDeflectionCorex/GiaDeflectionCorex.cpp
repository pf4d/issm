/*!\file GiaDeflectionCorex
 * \brief: GIA solution from Erik Ivins. 
 * Compute deflection wi from a single disk of radius re, load history hes for 
 * numtimes time steps. 
 */

#include "./GiaDeflectionCorex.h"

#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../InputUpdateFromConstantx/InputUpdateFromConstantx.h"

/*External blocks: {{{*/
struct blockp{
	double pset[7];
};

struct blocko{
	double rhoi;
};

struct blockrad{
	double distrad; 
};

struct blocks{
	double aswokm_w; 
	double aswokm_dwdt; 
};

extern "C" { 
	int distme_(int* pNtime,int* pNtimp,int* pNtimm,double* time,double* bi,double* dmi,double* zhload);

	int what0_(int* piedge,int* pNtimp,int* pNtimm,double* time,double* bi,double* dmi);
	extern struct blockp blockp_;
	extern struct blocko blocko_;
	extern struct blockrad blockrad_;
	extern struct blocks blocks_;
}

/*}}}*/

void GiaDeflectionCorex( IssmDouble* pwi, IssmDouble* pdwidt, GiaDeflectionCoreArgs* arguments){

	/*intermediary: */
	int i;

	/*output: */
	IssmDouble wi=0;
	IssmDouble dwidt=0;

	/*inputs: {{{*/
	/*constant: */
	int idisk=1; // disk #
	IssmDouble yts;

	/*coming from the model structure, runtime configurable:*/
	/*gia solution parameters: */
	int iedge=0; 

	/*gia inputs: */
	IssmDouble ri; //radial distance from center of disk to vertex  i
	IssmDouble re; //radius of disk
	IssmDouble* hes; //loading history (in ice thickness)
	IssmDouble* times; //loading history times
	int numtimes; //loading history length
	IssmDouble currenttime;
	int Ntime; // number of times with load history 
	int Ntimm; // Ntime-1 : for slope/y-cept of load segments 
	int Ntimp; // Ntime+1 : for evaluation time  
	IssmDouble* blockt_time=NULL;
	IssmDouble* blockt_bi=NULL;
	IssmDouble* blockt_dmi=NULL;
	IssmDouble* blocky_zhload=NULL;

	/*gia material parameters: */
	IssmDouble lithosphere_shear_modulus;
	IssmDouble lithosphere_density;
	IssmDouble mantle_shear_modulus;
	IssmDouble mantle_viscosity;
	IssmDouble mantle_density;
	IssmDouble lithosphere_thickness;

	/*ice properties: */
	IssmDouble rho_ice;

	/*some debug info: */
	int disk_id;

/*}}}*/

	/*Recover material parameters and loading history: see GiaDeflectionCoreArgs for more details {{{*/
	ri                        = arguments->ri;
	re                        = arguments->re;
	hes                       = arguments->hes;
	times                     = arguments->times;
	numtimes                  = arguments->numtimes;
	currenttime               = arguments->currenttime;
	lithosphere_shear_modulus = arguments->lithosphere_shear_modulus;
	lithosphere_density       = arguments->lithosphere_density;
	mantle_shear_modulus      = arguments->mantle_shear_modulus;
	mantle_viscosity          = arguments->mantle_viscosity;
	mantle_density            = arguments->mantle_density;
	lithosphere_thickness     = arguments->lithosphere_thickness;
	rho_ice                   = arguments->rho_ice;
	disk_id                   = arguments->idisk;
	iedge                     = arguments->iedge;
	yts                       = arguments->yts;

	/*}}}*/

	/*Modify inputs to match naruse code: */
	Ntime=numtimes;
	Ntimm=Ntime-1;
	Ntimp=Ntime+1;

	/*Prepare block inputs for fortran distme and what0 routines of the naruse code: {{{*/
	/*Now, let's set pset from the data that we got in input to GiaDeflectionCorex: */
	blockp_.pset[0]=lithosphere_thickness;
	blockp_.pset[1]=mantle_viscosity;
	blockp_.pset[2]=lithosphere_shear_modulus;
	blockp_.pset[3]=mantle_shear_modulus;
	blockp_.pset[4]=lithosphere_density;
	blockp_.pset[5]=mantle_density;
	blockp_.pset[6]=re;
	blocko_.rhoi=rho_ice; 

	/*loading history: */
	blocky_zhload=xNew<IssmDouble>(Ntime);
	for(i=0;i<Ntime;i++){
	blocky_zhload[i]=hes[i];
	}

	/*times in kyr: */
	blockt_time=xNew<IssmDouble>(Ntimp);
	for (i=0;i<Ntimp;i++){
		blockt_time[i]=times[i]/1000.0/yts; 
		if(i==numtimes-1)blockt_time[i]=times[numtimes-1]/1000.0/yts; // final loading time, same as evaluation time
		if(i==numtimes)blockt_time[i]=times[numtimes-1]/1000.0/yts;   // evaluation time
	}

	/*bi: */
	blockt_bi=xNew<IssmDouble>(Ntimm);

	/*dmi: */
	blockt_dmi=xNew<IssmDouble>(Ntimm);

	/*radial distance of i-th element: */
	blockrad_.distrad=ri/1000.0; // in km
	/*}}}*/

	/*Call distme driver: */
	distme_(&Ntime,&Ntimp,&Ntimm,blockt_time,blockt_bi,blockt_dmi,blocky_zhload); 

	/*Call what0 driver: */
	what0_(&iedge,&Ntimp,&Ntimm,blockt_time,blockt_bi,blockt_dmi); 

	/*output solution: */
	wi = blocks_.aswokm_w;
	dwidt = blocks_.aswokm_dwdt;
	*pwi=wi;
	*pdwidt=dwidt;

	/*Free ressources: */
	xDelete<IssmDouble>(blockt_time);
	xDelete<IssmDouble>(blockt_bi);
	xDelete<IssmDouble>(blockt_dmi);
	xDelete<IssmDouble>(blocky_zhload);

}
