/*!\file GEMB module from Alex Gardner.
 * \brief: calculates SMB 
 */

#include "./SurfaceMassBalancex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

const double Pi = 3.141592653589793;

void Gembx(FemModel* femmodel){  /*{{{*/

	for(int i=0;i<femmodel->elements->Size();i++){
        Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
        element->SmbGemb();
	}

} /*}}}*/
void GembgridInitialize(IssmDouble** pdz, int* psize, IssmDouble zTop, IssmDouble dzTop, IssmDouble zMax, IssmDouble zY){ /*{{{*/

	/* This file sets up the initial grid spacing and total grid depth.  The
	grid structure is set as constant grid length 'dzTop' for the top
	'zTop' meters of the model grid. Bellow 'zTop' the gid length increases
	linearly with depth */

	/*intermediary:*/
	IssmDouble dgpTop;
	int gpTop, gpBottom;
	int i;
	IssmDouble gp0,z0;
	IssmDouble* dzT=NULL;
	IssmDouble* dzB=NULL;

	/*output: */
	int m;
	IssmDouble* dz=NULL;

	//----------------------Calculate Grid Lengths------------------------------
	//calculate number of top grid points
	dgpTop = zTop/dzTop;

	//check to see if the top grid cell structure length (dzTop) goes evenly 
	//into specified top structure depth (zTop). Also make sure top grid cell
	//structure length (dzTop) is greater than 5 cm
	#ifndef _HAVE_ADOLC_  //avoid the round operation check!
	if (dgpTop != round(dgpTop)){ 
		_error_("top grid cell structure length does not go evenly into specified top structure depth, adjust dzTop or zTop");
	}
	#endif
	if(dzTop < 0.05){
		_printf_("initial top grid cell length (dzTop) is < 0.05 m");
	}
	gpTop=reCast<int,IssmDouble>(dgpTop);

	//initialize top grid depth vector
	dzT = xNew<IssmDouble>(gpTop); 
	for (i=0;i<gpTop;i++)dzT[i]=dzTop;

	//bottom grid cell depth = x*zY^(cells from to structure)
	//figure out the number of grid points in the bottom vector (not known a priori)
	gp0 = dzTop;
	z0 = zTop;
	gpBottom = 0;
	while (zMax > z0){
		gp0= gp0 * zY;
		z0 = z0 + gp0;
		gpBottom++;
	}
	//initialize bottom vectors
	dzB = xNewZeroInit<IssmDouble>(gpBottom);
	gp0 = dzTop;
	z0 = zTop;
	for(i=0;i<gpBottom;i++){
		gp0=gp0*zY;
		dzB[i]=gp0;
	}

	//combine top and bottom dz vectors
	dz = xNew<IssmDouble>(gpTop+gpBottom);
	for(i=0;i<gpTop;i++){
		dz[i]=dzT[i];
	}
	for(i=0;i<gpBottom;i++){
		dz[gpTop+i]=dzB[i];
	}

	/*Free ressouces:*/
	xDelete(dzT);
	xDelete(dzB);
	
	//---------NEED TO IMPLEMENT A PROPER GRID STRECHING ALGORITHM------------

	/*assign ouput pointers: */
	*pdz=dz; 
	*psize=gpTop+gpBottom;
} /*}}}*/ 
IssmDouble Marbouty(IssmDouble T, IssmDouble d, IssmDouble dT){ /*{{{*/

	// calculates grain growth according to Fig. 9 of Marbouty, 1980
	// ------NO GRAIN GROWTH FOR d > 400 kg m-3 because H is set to zero------
	// ---------------this is a major limitation of the model-------------------

	// initialize
	IssmDouble F = 0, H=0, G=0;
	const IssmDouble E = 0.09;        //[mm d-1] model time growth constant E
	// convert T from K to ºC
	T = T - 273.15;

	// temperature coefficient F
	if(T>-6.0) F =  0.7 + ((T/-6.0) * 0.3);
	if(T<=-6.0 && T>-22.0) F =  1 - ((T+6.0)/-16.0 * 0.8);
	if(T<=-22.0 && T>-40.0) F =  0.2 - ((T+22.0)/-18.0 * 0.2);

	// density coefficient F
	if(d<150.0) H=1.0;

	if(d>=150.0 & d<400.0) H = 1 - ((d-150.0)/250.0);

	// temperature gradient coefficient G
	if(dT >= 0.16 && dT < 0.25) G = ((dT - 0.16)/0.09) * 0.1;
	if(dT >= 0.25 && dT < 0.4)  G = 0.1 + (((dT - 0.25)/0.15) * 0.57);
	if(dT >= 0.4  && dT < 0.5)  G = 0.67 + (((dT - 0.4)/0.1) * 0.23);
	if(dT >= 0.5  && dT < 0.7)  G = 0.9 + (((dT - 0.5)/0.2) * 0.1);
	if(dT >= .7              )  G = 1.0;

	// grouped coefficient Q
	return F*H*G*E;

} /*}}}*/
void grainGrowth(IssmDouble* re, IssmDouble* gdn, IssmDouble* gsp, IssmDouble* T,IssmDouble* dz,IssmDouble* d, IssmDouble* W,IssmDouble smb_dt,int m,int aIdx,int sid){ /*{{{*/

	/*Created by: Alex S. Gardner, University of Alberta

	 *Description*: models the effective snow grain size

	 *Reference:*
	 DENDRITIC SNOW METAMORPHISM:
	 Brun, E., P. David, M. Sudul, and G. Brunot, 1992: A numerical model to
	 simulate snow-cover stratigraphy for operational avalanche forecasting.
	 Journal of Glaciology, 38, 13-22.

	 NONDENDRITIC SNOW METAMORPHISM:
	 Dry snow metamorphism:
	 Marbouty, D., 1980: An experimental study of temperature-gradient
	 metamorphism. Journal of Glaciology, 26, 303-312.

	 WET SNOW METAMORPHISM:
	 Brun, E., 1989: Investigation on wet-snow metamorphism in respect of
	 liquid-water content. Annals of Glaciology, 13, 22-26.

	 INPUTS
	 * T: grid cell temperature [k]
	 * dz: grid cell depth [m]
	 * d: grid cell density [kg m-3]
	 * W: water content [kg/m^2]
	 * re: effective grain size [mm]
	 * gdn: grain dentricity
	 * gsp: grain sphericity
	 * dt: time step of input data [s]

	 OUTPUTS
	 * re: effective grain size [mm]
	 * gdn: grain dentricity
	 * gsp: grain sphericity*/

	/*intermediary: */
	IssmDouble  dt;
	IssmDouble* gsz=NULL;
	IssmDouble* dT=NULL;
	IssmDouble* zGPC=NULL;
	IssmDouble* lwc=NULL;
	IssmDouble  Q=0;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   grain growth module\n");
	
	/*only when aIdx = 1 or 2 do we run grainGrowth: */
	if(aIdx!=1 & aIdx!=2){
		/*come out as we came in: */
		return;
	}

	/*Figure out grain size from effective grain radius: */
	gsz=xNew<IssmDouble>(m); for(int i=0;i<m;i++)gsz[i]=re[i]*2.0;

	/*Convert dt from seconds to day: */
	dt=smb_dt/86400.0;

	/*Determine liquid-water content in percentage: */
	lwc=xNew<IssmDouble>(m); for(int i=0;i<m;i++)lwc[i]= W[i] / (d[i]*dz[i])*100.0;

	//set maximum water content by mass to 9 percent (Brun, 1980)
	for(int i=0;i<m;i++)if(lwc[i]>9.0) lwc[i]=9.0;


	/* Calculate temperature gradiant across grid cells 
	 * Returns the avereage gradient across the upper and lower grid cell */

	//initialize
	dT=xNewZeroInit<IssmDouble>(m); 

	//depth of grid point center from surface
	zGPC=xNewZeroInit<IssmDouble>(m);  
	for(int i=0;i<m;i++){
		for (int j=0;j<=i;j++) zGPC[i]+=dz[j];
		zGPC[i]-=dz[i]/2;
	}

	// Take forward differences on left and right edges
	dT[0] = (T[1] - T[0])/(zGPC[1]-zGPC[0]);
	dT[m-1] = (T[m-1] - T[m-2])/(zGPC[m-1]-zGPC[m-2]);

	//Take centered differences on interior points
	for(int i=1;i<m-1;i++) dT[i] = (T[i+1]-T[i-1])/(zGPC[i+1]-zGPC[i-1]);

	// take absolute value of temperature gradient
	for(int i=0;i<m;i++)dT[i]=fabs(dT[i]);
	
	/*Snow metamorphism. Depends on value of snow dendricity and wetness of the snowpack: */
	for(int i=0;i<m;i++){
		if (gdn[i]>0){

			//_printf_("Dendritic dry snow metamorphism\n");

			//index for dentricity > 0 and T gradients < 5 degC m-1 and >= 5 degC m-1
			if(dT[i]<5.0){
				//determine coefficients
				IssmDouble A = - 2e8 * exp(-6e3 / T[i]) * dt;
				IssmDouble B = 1e9 * exp(-6e3 / T[i]) * dt;
				//new dentricity and sphericity for dT < 5 degC m-1
				gdn[i] += A;
				gsp[i] += B;
			}
			else{
				// new dendricity and sphericity for dT >= 5 degC m-1

				//determine coefficients
				IssmDouble C = (-2e8 * exp(-6e3 / T[i]) * dt) * pow(dT[i],.4);
				gdn[i] +=C;
				gsp[i] +=C;
			}

			// wet snow metamorphism
			if(W[i]>0.0){

				//_printf_("D}ritic wet snow metamorphism\n");

				//determine coefficient
				IssmDouble D = (1.0/16.0) * pow(lwc[i],3.0) * dt;

				// new dendricity and sphericity for wet snow
				gdn[i] -= D;
				gsp[i] += D;
			}

			// dendricity and sphericity can not be > 1 or < 0
			if (gdn[i]<0.0)gdn[i]=0.0;
			if (gsp[i]<0.0)gsp[i]=0.0;
			if (gdn[i]>1.0)gdn[i]=1.0;
			if (gsp[i]>1.0)gsp[i]=1.0;

			// determine new grain size (mm)
			gsz[i] = 0.1 + (1-gdn[i])*0.25 + (0.5-gsp[i])*0.1;

		}
		else{

			/*Dry snow metamorphism (Marbouty, 1980) grouped model coefficinets
			 *from Marbouty, 1980: Figure 9*/

			//_printf_("Nond}ritic snow metamorphism\n");
			Q = Marbouty(T[i],d[i],dT[i]);

			// calculate grain growth
			gsz[i] += Q* dt;

			//Wet snow metamorphism (Brun, 1989)
			if(W[i]>0.0){
				//_printf_("Nond}ritic wet snow metamorphism\n");
				//wet rate of change coefficient
				IssmDouble E = 1.28e-8 + (4.22e-10 * pow(lwc[i],3.0))* (dt *86400.0);   // [mm^3 s^-1]

				// calculate change in grain volume and convert to grain size
				gsz[i] = 2.0 * pow(3.0/(Pi * 4.0)*((4.0/ 3.0)*Pi*pow(gsz[i]/2.0,3.0) + E),1.0/3.0);

			}

			// grains with sphericity == 1 can not have grain sizes > 2 mm (Brun, 1992)
			if (gsp[i]==1.0 && gsz[i]>2.0) gsz[i]=2.0;

			// grains with sphericity == 0 can not have grain sizes > 5 mm (Brun, 1992)
			if (gsp[i]!=1.0 && gsz[i]>5.0) gsz[i]=5.0;
		}

		//convert grain size back to effective grain radius: 
		re[i]=gsz[i]/2.0;
	}
	
	/*Free ressources:*/
	xDelete<IssmDouble>(gsz);
	xDelete<IssmDouble>(dT);
	xDelete<IssmDouble>(zGPC);
	xDelete<IssmDouble>(lwc);

}  /*}}}*/
void albedo(IssmDouble* a, int aIdx, IssmDouble* re, IssmDouble* d, IssmDouble cldFrac, IssmDouble aIce, IssmDouble aSnow, IssmDouble* TK, IssmDouble* W, IssmDouble P, IssmDouble EC, IssmDouble t0wet, IssmDouble t0dry, IssmDouble K, IssmDouble dt, int m,int sid) { /*{{{*/

	//// Calculates Snow, firn and ice albedo as a function of:
	//   1 : effective grain radius (Gardner & Sharp, 2009)
	//   2 : effective grain radius (Brun et al., 2009)
	//   3 : density and cloud amount (Greuell & Konzelmann, 1994)
	//   4 : exponential time decay & wetness (Bougamont & Bamber, 2005)

	//// Inputs
	// aIdx      = albedo method to use

	// Methods 1 & 2
	//   re      = surface effective grain radius [mm]

	// Methods 3
	//   d       = snow surface density [kg m-3]
	//   n       = cloud amount
	//   aIce    = albedo of ice
	//   aSnow   = albedo of fresh snow

	// Methods 4
	//   aIce    = albedo of ice
	//   aSnow   = albedo of fresh snow
	//   a       = grid cell albedo from prevous time step;
	//   T       = grid cell temperature [k]
	//   W       = pore water [kg]
	//   P       = precipitation [mm w.e.] or [kg m-3]
	//   EC      = surface evaporation (-) condensation (+) [kg m-2]
	//   t0wet   = time scale for wet snow (15-21.9) [d]
	//   t0dry   = warm snow timescale [15] [d]
	//   K       = time scale temperature coef. (7) [d]
	//   dt      = time step of input data [s]

	//// Output
	//   a       = grid cell albedo 

	//// Usage 
	// Method 1
	// a = albedo(1, 0.1); 

	// Method 4
	// a = albedo(4, [], [], [], 0.48, 0.85, [0.8 0.5 ... 0.48], ...
	//   [273 272.5 ... 265], [0 0.001 ... 0], 0, 0.01, 15, 15, 7, 3600)

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   albedo module\n");
        
	//some constants:
	const IssmDouble dSnow = 300;   // density of fresh snow [kg m-3]       
	const IssmDouble dIce = 910;    // density of ice [kg m-3]

	if(aIdx==1){ 
        //function of effective grain radius
        
        //convert effective radius to specific surface area [cm2 g-1]
        IssmDouble S = 3.0 / (.091 * re[0]);
        
        //determine broadband albedo
        a[0]= 1.48 - pow(S,-.07);
	}
	else if(aIdx==2){
		
        // Spectral fractions  (Lefebre et al., 2003)
        // [0.3-0.8um 0.8-1.5um 1.5-2.8um]
        
        IssmDouble sF[3] = {0.606, 0.301, 0.093};
        
        // convert effective radius to grain size in meters
        IssmDouble gsz = (re[0] * 2) / 1000.0;
        
        // spectral range:
        // 0.3 - 0.8um
        IssmDouble a0 = fmin(0.98, 1 - 1.58 *pow(gsz,0.5));
        // 0.8 - 1.5um
        IssmDouble a1 = fmax(0, 0.95 - 15.4 *pow(gsz,0.5));
        // 1.5 - 2.8um
        IssmDouble a2 = fmax(0.127, 0.88 + 346.3*gsz - 32.31*pow(gsz,0.5));
        
        // broadband surface albedo
        a[0] = sF[0]*a0 + sF[1]*a1 + sF[2]*a2;

	}
	else if(aIdx==3){
		
        // a as a function of density
        
        // calculate albedo
        a[0] = aIce + (d[0] - dIce)*(aSnow - aIce) / (dSnow - dIce) + (0.05 * (cldFrac - 0.5));
	}
	else if(aIdx==4){
		
        // exponential time decay & wetness
        
        // change in albedo with time:
        //   (d_a) = (a - a_old)/(t0)
        // where: t0 = timescale for albedo decay
        
        dt = dt / 86400;    // convert from [s] to [d]
        
        // initialize variables
        IssmDouble* t0=xNew<IssmDouble>(m);
        IssmDouble* T=xNew<IssmDouble>(m);
        IssmDouble* t0warm=xNew<IssmDouble>(m);
        IssmDouble* d_a=xNew<IssmDouble>(m);
        
        // specify constants
        // a_wet = 0.15;        // water albedo (0.15)
        // a_new = aSnow        // new snow albedo (0.64 - 0.89)
        // a_old = aIce;        // old snow/ice albedo (0.27-0.53)
        // t0_wet = t0wet;      // time scale for wet snow (15-21.9) [d]
        // t0_dry = t0dry;      // warm snow timescale [15] [d]
        // K = 7                // time scale temperature coef. (7) [d]
        // W0 = 300;            // 200 - 600 [mm]
        const IssmDouble z_snow = 15;            // 16 - 32 [mm]
        
        // determine timescale for albedo decay
        for(int i=0;i<m;i++)if(W[i]>0)t0[i]=t0wet; // wet snow timescale
        for(int i=0;i<m;i++)T[i]=TK[i] - 273.15; // change T from K to degC
        for(int i=0;i<m;i++) t0warm[i]= fabs(T[i]) * K + t0dry; //// 'warm' snow timescale
        for(int i=0;i<m;i++)if(W[i]==0.0 && T[i]>=-10)t0[i]= t0warm[i];
        for(int i=0;i<m;i++)if(T[i]<-10) t0[i] =  10 * K + t0dry; // 'cold' snow timescale
        
        // calculate new albedo
        for(int i=0;i<m;i++)d_a[i] = (a[i] - aIce) / t0[i] * dt;           // change in albedo
        for(int i=0;i<m;i++)a[i] -= d_a[i];                            // new albedo
        
        // modification of albedo due to thin layer of snow or solid
        // condensation (deposition) at the surface surface
        
        // check if condensation occurs & if it is deposited in solid phase
        if ( EC > 0 && T[0] < 0) P = P + (EC/dSnow) * 1000;  // add cond to precip [mm]
        
        a[0] = aSnow - (aSnow - a[0]) * exp(-P/z_snow);
        
        //----------THIS NEEDS TO BE IMPLEMENTED AT A LATER DATE------------
        // modification of albedo due to thin layer of water on the surface
        // a_surf = a_wet - (a_wet - a_surf) * exp(-W_surf/W0);

        /*Free ressources:*/
        xDelete<IssmDouble>(t0);
        xDelete<IssmDouble>(T);
        xDelete<IssmDouble>(t0warm);
        xDelete<IssmDouble>(d_a);

	}
	else _error_("albedo method switch should range from 1 to 4!");
	
	// Check for erroneous values
	if (a[0] > 1) _printf_("albedo > 1.0\n");
	else if (a[0] < 0) _printf_("albedo is negative\n");
	else if (xIsNan(a[0])) _error_("albedo == NAN\n");
}  /*}}}*/
void thermo(IssmDouble* pEC, IssmDouble* T, IssmDouble* dz, IssmDouble* d, IssmDouble* swf, IssmDouble dlwrf, IssmDouble Ta, IssmDouble V, IssmDouble eAir, IssmDouble pAir, IssmDouble Ws, IssmDouble dt0, int m, IssmDouble Vz, IssmDouble Tz,int sid) { /*{{{*/

	/* ENGLACIAL THERMODYNAMICS*/
	 
	/* Description: 
	   computes new temperature profile accounting for energy absorption and 
	   thermal diffusion.*/

	// INPUTS
	//  T: grid cell temperature [k]
	//  dz: grid cell depth [m]
	//  d: grid cell density [kg m-3]
	//  swf: shortwave radiation fluxes [W m-2]
	//  dlwrf: downward longwave radiation fluxes [W m-2]
	//  Ta: 2 m air temperature
	//  V:  wind velocity [m s-1]
	//  eAir: screen level vapor pressure [Pa]
	//  Ws: surface water content [kg]
	//  dt0: time step of input data [s]
	//  elev: surface elevation [m a.s.l.] 
	//  Vz: air temperature height above surface [m]
	//  Tz: wind height above surface [m]

	// OUTPUTS
	// T: grid cell temperature [k]
	// EC: evaporation/condensation [kg]
	
	/*intermediary: */
	IssmDouble* K = NULL;
	IssmDouble* KU = NULL;
	IssmDouble* KD = NULL;
	IssmDouble* KP = NULL;
	IssmDouble* Au = NULL;
	IssmDouble* Ad = NULL;
	IssmDouble* Ap = NULL;
	IssmDouble* Nu = NULL;
	IssmDouble* Nd = NULL;
	IssmDouble* Np = NULL;
	IssmDouble* dzU = NULL;
	IssmDouble* dzD = NULL;
	IssmDouble* sw = NULL;
	IssmDouble* dT_sw = NULL;
	IssmDouble* lw = NULL;
	IssmDouble* T0 = NULL;
	IssmDouble* Tu = NULL;
	IssmDouble* Td = NULL;

	IssmDouble  z0;	
	IssmDouble  dt;
	IssmDouble max_fdt=0;
	IssmDouble  Ts=0;
	IssmDouble  L;
	IssmDouble  eS;
	IssmDouble  Ri=0;
	IssmDouble  coefM;
	IssmDouble  coefH;
	IssmDouble An;
	IssmDouble C;
	IssmDouble shf;
	IssmDouble SB;
	IssmDouble CI; 
	IssmDouble ds;
	IssmDouble dAir;
	IssmDouble TCs;
	IssmDouble lhf;
	IssmDouble EC_day;
	IssmDouble dT_turb;
	IssmDouble turb;
	IssmDouble ulw;
	IssmDouble dT_ulw;
	IssmDouble dlw;
	IssmDouble dT_dlw;
	
	/*outputs:*/
	IssmDouble EC;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   thermal module\n");

	// INITIALIZE
	CI = 2102;          // heat capacity of snow/ice (J kg-1 k-1)
	// CA = 1005;                  // heat capacity of air (J kg-1 k-1)
	// LF = 0.3345E6;              // latent heat of fusion(J kg-1)
	// LV = 2.495E6;               // latent heat of vaporization(J kg-1)
	// dIce = 910;                 // density of ice [kg m-3]
	// dSnow = 300;                // density of snow [kg m-3]
	SB = 5.67E-8;       // Stefan-Boltzmann constant [W m-2 K-4]

	ds = d[0];      // density of top grid cell

	// calculated air density [kg/m3]
	dAir = 0.029 * pAir /(8.314 * Ta);

	// thermal capacity of top grid cell [J/k]
	TCs = d[0]*dz[0]*CI; 

	//initialize Evaporation - Condenstation 
	EC = 0;
	
	// check if all SW applied to surface or distributed throught subsurface
	// swIdx = length(swf) > 1

	// SURFACE ROUGHNESS (Bougamont, 2006)
	// wind/temperature surface roughness height [m]
	if (ds < 910 && Ws == 0) z0 = 0.00012;       // 0.12 mm for dry snow
	else if (ds >= 910) z0 = 0.0032;             // 3.2 mm for ice
	else z0 = 0.0013;                            // 1.3 mm for wet snow

	// if V = 0 goes to infinity therfore if V = 0 change
	if(V<.01)V=.01;
	
	// Bulk-transfer coefficient for turbulent fluxes
	An =  pow(0.4,2) / pow(log(Tz/z0),2);     // Bulk-transfer coefficient
	C = An * dAir * V;                        // shf & lhf common coefficient

	// THERMAL CONDUCTIVITY (Sturm, 1997: J. Glaciology)
	// calculate new K profile [W m-1 K-1]

	// initialize conductivity
	K= xNewZeroInit<IssmDouble>(m);

	// for snow and firn (density < 910 kg m-3) (Sturn et al, 1997)
	for(int i=0;i<m;i++) if(d[i]<910) K[i] = 0.138 - 1.01E-3 * d[i] + 3.233E-6 * (pow(d[i],2));

	// for ice (density >= 910 kg m-3)
	for(int i=0;i<m;i++) if(d[i]>=910) K[i] = 9.828 * exp(-5.7E-3*T[i]);

	// THERMAL DIFFUSION COEFFICIENTS
 
	// A discretization scheme which truncates the Taylor-Series expansion
	// after the 3rd term is used. See Patankar 1980, Ch. 3&4
 
	// discretized heat equation:
 
	//                 Tp = (Au*Tu° + Ad*Td° + (Ap-Au-Ad)Tp° + S) / Ap
 
	// where neighbor coefficients Au, Ap, & Ad are
 
	//                   Au = [dz_u/2KP + dz_p/2KE]^-1
	//                   Ad = [dz_d/2KP + dz_d/2KD]^-1 
	//                   Ap = d*CI*dz/Dt 

	// and u & d represent grid points up and down from the center grid point 
	// p and // u & d represent grid points up and down from the center grid 
	// point p and ° identifies previous time step values. S is a source term.

	// u, d, and p conductivities
	KU = xNew<IssmDouble>(m);
	KD = xNew<IssmDouble>(m);
	KP = xNew<IssmDouble>(m);
	
	KU[0] = UNDEF;
	KD[m-1] = UNDEF;
	for(int i=1;i<m;i++) KU[i]= K[i-1];
	for(int i=0;i<m-1;i++) KD[i] = K[i+1];
	for(int i=0;i<m;i++) KP[i] = K[i];

	// determine u, d & p cell widths
	dzU = xNew<IssmDouble>(m);
	dzD = xNew<IssmDouble>(m);
	dzU[0]=UNDEF;
	dzD[m-1]=UNDEF;
	
	for(int i=1;i<m;i++) dzU[i]= dz[i-1];
	for(int i=0;i<m-1;i++) dzD[i] = dz[i+1];

	// determine minimum acceptable delta t (diffusion number > 1/2) [s]
	dt=1e12; 
	for(int i=0;i<m;i++)dt = fmin(dt,CI * pow(dz[i],2) * d[i]  / (3 * K[i]));

	// smallest possible even integer of 60 min where diffusion number > 1/2
	// must go evenly into one hour or the data frequency if it is smaller

	// all integer factors of the number of second in a day (86400 [s])
	int f[45] = {1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 30, 36, 40, 45, 48, 50, 60,
    72, 75, 80, 90, 100, 120, 144, 150, 180, 200, 225, 240, 300, 360, 400, 450, 600, 720, 900, 1200, 1800, 3600};

	// return the min integer factor that is < dt
	max_fdt=f[0];
	for(int i=0;i<45;i++){
		if (f[i]<dt)if(f[i]>=max_fdt)max_fdt=f[i];
	}
	dt=max_fdt;
	
	// determine mean (harmonic mean) of K/dz for u, d, & p
	Au = xNew<IssmDouble>(m);
	Ad = xNew<IssmDouble>(m);
	Ap = xNew<IssmDouble>(m);
	for(int i=0;i<m;i++){
		Au[i] = pow((dzU[i]/2/KP[i] + dz[i]/2/KU[i]),-1);
		Ad[i] = pow((dzD[i]/2/KP[i] + dz[i]/2/KD[i]),-1);
		Ap[i] = (d[i]*dz[i]*CI)/dt;
	}
	
	// create "neighbor" coefficient matrix
	Nu = xNew<IssmDouble>(m);
	Nd = xNew<IssmDouble>(m);
	Np = xNew<IssmDouble>(m);
	for(int i=0;i<m;i++){
		Nu[i] = Au[i] / Ap[i];
		Nd[i] = Ad[i] / Ap[i];
		Np[i]= 1 - Nu[i] - Nd[i];
	}
	
	// specify boundary conditions: constant flux at bottom
	Nu[m-1] = 0;
	Np[m-1] = 1;
	
	// zero flux at surface
	Np[0] = 1 - Nd[0];
	
	// Create neighbor arrays for diffusion calculations instead of a tridiagonal matrix
	Nu[0] = 0;
	Nd[m-1] = 0;
	
	/* RADIATIVE FLUXES*/

	// energy supplied by shortwave radiation [J]
	sw = xNew<IssmDouble>(m);
	for(int i=0;i<m;i++) sw[i]= swf[i] * dt;
	
	// temperature change due to SW
	dT_sw = xNew<IssmDouble>(m);
	for(int i=0;i<m;i++) dT_sw[i]= sw[i] / (CI * d[i] * dz[i]);

	// Upward longwave radiation flux is calculated from the snow surface
	// temperature which is set equal to the average temperature of the
	// top grid cells.

	// energy supplied by downward longwave radiation to the top grid cell [J]
	dlw = dlwrf * dt;

	// temperature change due to dlw_surf
	dT_dlw = dlw / TCs;

	// PREALLOCATE ARRAYS BEFORE LOOP FOR IMPROVED PERFORMANCE
	T0 = xNewZeroInit<IssmDouble>(m+2);
	Tu=xNew<IssmDouble>(m);
	Td=xNew<IssmDouble>(m);

	/* CALCULATE ENERGY SOURCES AND DIFFUSION FOR EVERY TIME STEP [dt]*/
	for (IssmDouble i=1;i<=dt0;i+=dt){

		// PART OF ENERGY CONSERVATION CHECK
		// store initial temperature
		//T_init = T;
    
		// calculate temperature of snow surface (Ts)
		// when incoming SW radition is allowed to penetrate the surface,
		// the modeled energy balance becomes very sensitive to how Ts is
		// calculated.  The estimated enegy balance & melt are significanly
		// less when Ts is taken as the mean of the x top grid cells.
		Ts = (T[0] + T[1])/2;
		Ts = fmin(273.15,Ts);    // don't allow Ts to exceed 273.15 K (0°C)
		
		//TURBULENT HEAT FLUX
    
		// Monin–Obukhov Stability Correction
		// Reference:
		// Ohmura, A., 1982: Climate and Energy-Balance on the Arctic Tundra.
		// Journal of Climatology, 2, 65-84.

		// calculate the Bulk Richardson Number (Ri)
		Ri = (2*9.81* (Vz - z0) * (Ta - Ts)) / ((Ta + Ts)* pow(V,2.0));
		
		// calculate Monin-–Obukhov stability factors 'coefM' and 'coefH'
    
		// do not allow Ri to exceed 0.19
		Ri = fmin(Ri, 0.19);

		// calculate momentum 'coefM' stability factor
		if (Ri > 0){
			// if stable
			coefM = 1/(1-5.2*Ri);
		}
		else {
			coefM =pow (1-18*Ri,-0.25);
		}
		
		// calculate heat/wind 'coef_H' stability factor
		if (Ri < -0.03) coefH = 1.3 * coefM;
		else coefH = coefM;
		
		//// Sensible Heat
		// calculate the sensible heat flux [W m-2](Patterson, 1998)
		shf = C * 1005 * (Ta - Ts);

		// adjust using Monin–Obukhov stability theory
		shf = shf / (coefM * coefH);

		//// Latent Heat
		// determine if snow pack is melting & calcualte surface vapour pressure over ice or liquid water
		if (Ts >= 273.15){
			L = 2.495E6;

			// for an ice surface Murphy and Koop, 2005 [Equation 7]
			eS = exp(9.550426 - 5723.265/Ts + 3.53068 * log(Ts) - 0.00728332 * Ts);
		}
		else{
			L = 2.8295E6; // latent heat of sublimation for liquid surface (assume liquid on surface when Ts == 0 deg C)
			// Wright (1997), US Meteorological Handbook from Murphy and Koop, 2005 Appendix A
			eS = 611.21 * exp(17.502 * (Ts - 273.15) / (240.97 + Ts - 273.15));
		}
		
		// Latent heat flux [W m-2]
		lhf = C * L * (eAir - eS) * 0.622 / pAir;
		
		// adjust using Monin–Obukhov stability theory (if lhf '+' then there is energy and mass gained at the surface, 
		// if '-' then there is mass and energy loss at the surface.
		lhf = lhf / (coefM * coefH);

		//mass loss (-)/acreation(+) due to evaporation/condensation [kg]
		EC_day = lhf * 86400 / L;

		// temperature change due turbulent fluxes
		turb = (shf + lhf)* dt;
		dT_turb = turb  / TCs;

		// upward longwave contribution
		ulw = - SB * pow(Ts,4.0) * dt;
		dT_ulw = ulw / TCs;
		
		// new grid point temperature
    
		//SW penetrates surface
		for(int j=0;j<m;j++) T[j] = T[j] + dT_sw[j];
		T[0] = T[0] + dT_dlw + dT_ulw + dT_turb;
		
		// temperature diffusion
		for(int j=0;j<m;j++)T0[1+j]=T[j];
		for(int j=0;j<m;j++) Tu[j] = T0[j];
		for(int j=0;j<m;j++) Td[j] = T0[2+j];
		for(int j=0;j<m;j++) T[j] = (Np[j] * T[j]) + (Nu[j] * Tu[j]) + (Nd[j] * Td[j]);
		
		// calculate cumulative evaporation (+)/condensation(-)
		EC = EC + (EC_day/86400)*dt;
    
		/* CHECK FOR ENERGY (E) CONSERVATION [UNITS: J]
		//energy flux across lower boundary (energy supplied by underling ice)
		base_flux = Ad(-1)*(T_init()-T_init(-1)) * dt;

		E_used = sum((T - T_init) * (d*dz*CI));
		E_sup = ((sum(swf)  * dt) + dlw + ulw + turb + base_flux);

		E_diff = E_used - E_sup;

		if abs(E_diff) > 1E-6 || isnan(E_diff)
		disp(T(1))
		_error_("energy not conserved in thermodynamics equations");
		*/
	}
	
	/*Free ressources:*/
	xDelete<IssmDouble>(K);
	xDelete<IssmDouble>(KU);
	xDelete<IssmDouble>(KD);
	xDelete<IssmDouble>(KP);
	xDelete<IssmDouble>(Au);
	xDelete<IssmDouble>(Ad);
	xDelete<IssmDouble>(Ap);
	xDelete<IssmDouble>(Nu);
	xDelete<IssmDouble>(Nd);
	xDelete<IssmDouble>(Np);
	xDelete<IssmDouble>(dzU);
	xDelete<IssmDouble>(dzD);
	xDelete<IssmDouble>(sw);
	xDelete<IssmDouble>(dT_sw);
	xDelete<IssmDouble>(lw);
	xDelete<IssmDouble>(T0);
	xDelete<IssmDouble>(Tu);
	xDelete<IssmDouble>(Td);


	/*Assign output pointers:*/
	*pEC=EC;

}  /*}}}*/
void shortwave(IssmDouble** pswf, int swIdx, int aIdx, IssmDouble dsw, IssmDouble as, IssmDouble* d, IssmDouble* dz, IssmDouble* re, int m, int sid){ /*{{{*/

	// DISTRIBUTES ABSORBED SHORTWAVE RADIATION WITHIN SNOW/ICE

	// swIdx = 0 : all absorbed SW energy is assigned to the top grid cell

	// swIdx = 1 : absorbed SW is distributed with depth as a function of:
	//   1 : snow density (taken from Bassford, 2004)
	//   2 : grain size in 3 spectral bands (Brun et al., 1992)

	// Inputs
	//   swIdx   = shortwave allowed to penetrate surface (0 = No, 1 = Yes)
	//   aIdx    = method for calculating albedo (1-4)
	//   dsw     = downward shortwave radiative flux [w m-2]
	//   as      = surface albedo
	//   d       = grid cell density [kg m-3]
	//   dz      = grid cell depth [m]
	//   re      = grid cell effective grain radius [mm]

	// Outputs
	//   swf     = absorbed shortwave radiation [W m-2]
	//
	
	/*outputs: */
	IssmDouble* swf=NULL;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   shortwave module\n");

	/*Initialize and allocate: */
	swf=xNewZeroInit<IssmDouble>(m);


	// SHORTWAVE FUNCTION
	if (swIdx == 0) {// all sw radation is absorbed in by the top grid cell
        
		// calculate surface shortwave radiation fluxes [W m-2]
		swf[0] = (1 - as) * dsw;
	}
	else{ // sw radation is absorbed at depth within the glacier

		if (aIdx == 2){    // function of effective radius (3 spectral bands)

			IssmDouble * gsz=NULL;
			IssmDouble * B1_cum=NULL;
			IssmDouble * B2_cum=NULL;
			IssmDouble* h =NULL;
			IssmDouble* B1 =NULL;
			IssmDouble* B2 =NULL;
			IssmDouble* exp1 = NULL;
			IssmDouble* exp2 = NULL;
			IssmDouble*  Qs1 = NULL;
			IssmDouble*  Qs2 = NULL;

			// convert effective radius [mm] to grain size [m]
			gsz=xNew<IssmDouble>(m);
			for(int i=0;i<m;i++) gsz[i]= (re[i] * 2) / 1000;

			// Spectral fractions [0.3-0.8um 0.8-1.5um 1.5-2.8um]
			// (Lefebre et al., 2003)
			IssmDouble sF[3] = {0.606, 0.301, 0.093};

			// initialize variables
			B1_cum=xNew<IssmDouble>(m+1);
			B2_cum=xNew<IssmDouble>(m+1);
			for(int i=0;i<m+1;i++){
				B1_cum[i]=1;
				B2_cum[i]=1;
			}


			// spectral albedos:
			// 0.3 - 0.8um
			IssmDouble a0 = fmin(0.98, 1 - 1.58 *pow(gsz[0],0.5));
			// 0.8 - 1.5um
			IssmDouble a1 = fmax(0, 0.95 - 15.4 *pow(gsz[0],0.5));
			// 1.5 - 2.8um
			IssmDouble a2 = fmax(0.127, 0.88 + 346.3*gsz[0] - 32.31*pow(gsz[0],0.5));

			// separate net shortwave radiative flux into spectral ranges
			IssmDouble swfS[3];
			swfS[0] = (sF[0] * dsw) * (1 - a0);
			swfS[1] = (sF[1] * dsw) * (1 - a1);
			swfS[2] = (sF[2] * dsw) * (1 - a2);

			// absorption coeficient for spectral range:
			h =xNew<IssmDouble>(m);
			B1 =xNew<IssmDouble>(m);
			B2 =xNew<IssmDouble>(m);
			for(int i=0;i<m;i++) h[i]= d[i] /(pow(gsz[i],0.5));
			for(int i=0;i<m;i++) B1[i] = .0192 * h[i];                 // 0.3 - 0.8um
			for(int i=0;i<m;i++) B2[i]= .1098 * h[i];                 // 0.8 - 1.5um
			// B3 = +inf                     // 1.5 - 2.8um

			// cumulative extinction factors
			exp1 = xNew<IssmDouble>(m); 
			exp2 = xNew<IssmDouble>(m); 
			for(int i=0;i<m;i++) exp1[i]=exp(-B1[i]*dz[i]);
			for(int i=0;i<m;i++) exp2[i]=exp(-B2[i]*dz[i]);

			for(int i=0;i<m;i++){
				IssmDouble cum1=exp1[0];
				IssmDouble cum2=exp2[0];
				for(int j=1;j<=i;j++){
					cum1 = cum1*exp1[j];
					cum2 = cum2*exp2[j];
				}
				B1_cum[i+1]=cum1;
				B2_cum[i+1]=cum2;
			}


			// flux across grid cell boundaries
			Qs1 = xNew<IssmDouble>(m+1);
			Qs2 = xNew<IssmDouble>(m+1);
			for(int i=0;i<m+1;i++){
				Qs1[i] = swfS[0] * B1_cum[i];
				Qs2[i] = swfS[1] * B2_cum[i];
			}

			// net energy flux to each grid cell
			for(int i=0;i<m;i++) swf[i]= (Qs1[i]-Qs1[i+1]) + (Qs2[i]-Qs2[i+1]);

			// add flux absorbed at surface
			swf[0] = swf[0]+ swfS[2];

			/*Free ressources: */
			xDelete<IssmDouble>(gsz);
			xDelete<IssmDouble>(B1_cum);
			xDelete<IssmDouble>(B2_cum);
			xDelete<IssmDouble>(h);
			xDelete<IssmDouble>(B1);
			xDelete<IssmDouble>(B2);
			xDelete<IssmDouble>(exp1);
			xDelete<IssmDouble>(exp2);
			xDelete<IssmDouble>(Qs1);
			xDelete<IssmDouble>(Qs2);
			
			
		}
		else{  //function of grid cell density

			/*intermediary: */
			IssmDouble* B_cum = NULL;
			IssmDouble* exp_B = NULL;
			IssmDouble* Qs    = NULL;
			IssmDouble* B    = NULL;

			// fraction of sw radiation absorbed in top grid cell (wavelength > 0.8um)
			IssmDouble SWs = 0.36;

			// SWs and SWss coefficients need to be better constranted. Greuell
			// and Konzelmann 1994 used SWs = 0.36 and SWss = 0.64 as this the
			// the // of SW radiation with wavelengths > and < 800 nm
			// respectively.  This, however, may not account for the fact that
			// the albedo of wavelengths > 800 nm has a much lower albedo.

			// calculate surface shortwave radiation fluxes [W m-2]
			IssmDouble swf_s = SWs * (1 - as) * dsw;

			// calculate surface shortwave radiation fluxes [W m-2]
			IssmDouble swf_ss = (1-SWs) * (1 - as) * dsw;

			// SW allowed to penetrate into snowpack
			IssmDouble Bs = 10;    // snow SW extinction coefficient [m-1] (Bassford,2006)
			IssmDouble Bi = 1.3;   // ice SW extinction coefficient [m-1] (Bassford,2006)

			// calculate extinction coefficient B [m-1] vector
			B=xNew<IssmDouble>(m);
			for(int i=0;i<m;i++) B[i] = Bs + (300 - d[i]) * ((Bs - Bi)/(910 - 300));

			// cumulative extinction factor
			B_cum = xNew<IssmDouble>(m+1);
			exp_B = xNew<IssmDouble>(m);
			for(int i=0;i<m;i++)exp_B[i]=exp(-B[i]*dz[i]);

			B_cum[0]=1;
			for(int i=0;i<m;i++){
				IssmDouble cum_B=exp_B[0];
				for(int j=1;j<=i;j++) cum_B=cum_B*exp_B[j];
				B_cum[i+1]=  cum_B;
			}

			// flux across grid cell boundaries
			Qs=xNew<IssmDouble>(m+1);
			for(int i=0;i<m+1;i++) Qs[i] = swf_ss * B_cum[i];

			// net energy flux to each grid cell
			for(int i=0;i<m;i++) swf[i] = (Qs[i]-Qs[i+1]);

			// add flux absorbed at surface
			swf[0] += swf_s;

			/*Free ressources:*/
			xDelete<IssmDouble>(B_cum);
			xDelete<IssmDouble>(exp_B);
			xDelete<IssmDouble>(Qs);
			xDelete<IssmDouble>(B);
		}
	}
	/*Assign output pointers: */
	*pswf=swf;

} /*}}}*/ 
void accumulation(IssmDouble** pT, IssmDouble** pdz, IssmDouble** pd, IssmDouble** pW, IssmDouble** pa, IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, int* pm, IssmDouble T_air, IssmDouble P, IssmDouble dzMin, IssmDouble aSnow, int sid){ /*{{{*/

	// Adds precipitation and deposition to the model grid

	// Author: Alex Gardner, University of Alberta
	// Date last modified: JAN, 2008

	/* Description:
	   adjusts the properties of the top grid cell to account for accumulation
	   T_air & T = Air and top grid cell temperatures [K]
	   dz = topgrid cell length [m]
	   d = density of top grid gell [kg m-3]
	   P = precipitation [mm w.e.] or [kg m-3]
	   re = effective grain radius [mm]
	   gdn = grain dentricity
	   gsp = grain sphericity*/

	// MAIN FUNCTION
	// specify constants
	const IssmDouble dIce = 910;     // density of ice [kg m-3]
	const IssmDouble dSnow = 150;    // density of snow [kg m-3]
	const IssmDouble reNew = 0.1;    // new snow grain size [mm]
	const IssmDouble gdnNew = 1;     // new snow dendricity 
	const IssmDouble gspNew = 0.5;   // new snow sphericity 

	/*intermediary: */
	IssmDouble* mInit=NULL;
	bool        top=true;
	IssmDouble  mass, massinit, mass_diff;

	/*output: */
	IssmDouble* T=NULL;
	IssmDouble* dz=NULL;
	IssmDouble* d=NULL;
	IssmDouble* W=NULL;
	IssmDouble* a=NULL;
	IssmDouble* re=NULL;
	IssmDouble* gdn=NULL;
	IssmDouble* gsp=NULL;
	int         m;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   accumulation module\n");

	/*Recover pointers: */
	T=*pT;
	dz=*pdz;
	d=*pd;
	W=*pW;
	a=*pa;
	re=*pre;
	gdn=*pgdn;
	gsp=*pgsp;
	m=*pm;

	// determine initial mass
	mInit=xNew<IssmDouble>(m);
	for(int i=0;i<m;i++) mInit[i]= d[i] * dz[i];
	massinit=0; for(int i=0;i<m;i++)massinit+=mInit[i];

	if (P > 0){
			

		if (T_air <= 273.15){ // if snow

			IssmDouble  z_snow = P/dSnow;               // depth of snow

			// if snow depth is greater than specified min dz, new cell created
			if (z_snow > dzMin){

				newcell(&T,T_air,top,m); //new cell T
				newcell(&dz,z_snow,top,m); //new cell dz
				newcell(&d,dSnow,top,m); //new cell d
				newcell(&W,0,top,m); //new cell W
				newcell(&a,aSnow,top,m); //new cell a
				newcell(&re,reNew,top,m); //new cell grain size
				newcell(&gdn,gdnNew,top,m); //new cell grain dendricity
				newcell(&gsp,gspNew,top,m); //new cell grain sphericity
				m=m+1;
			}
			else { // if snow depth is less than specified minimum dz snow

				IssmDouble mass = mInit[0] + P;         // grid cell adjust mass

				dz[0] = dz[0] + P/dSnow;    // adjust grid cell depth      
				d[0] = mass / dz[0];    // adjust grid cell density

				// adjust variables as a linearly weighted function of mass
				// adjust temperature (assume P is same temp as air)
				T[0] = (T_air * P + T[0] * mInit[0])/mass;

				// adjust a, re, gdn & gsp
				a[0] = (aSnow * P + a[0] * mInit[0])/mass;
				re[0] = (reNew * P + re[0] * mInit[0])/mass;
				gdn[0] = (gdnNew * P + gdn[0] * mInit[0])/mass;
				gsp[0] = (gspNew * P + gsp[0] * mInit[0])/mass;
			}
		}
		else{ // if rain    

			/*rain is added by increasing the mass and temperature of the ice
			  of the top grid cell.  Temperatures are set artifically high to
			  account for the latent heat of fusion.  This is the same as
			  directly adding liquid water to the the snow pack surface but
			  makes the numerics easier.*/

			IssmDouble LF = 0.3345E6;  // latent heat of fusion(J kg-1)
			IssmDouble CI = 2102;      // specific heat capacity of snow/ice (J kg-1 k-1)

			// grid cell adjust mass
			mass = mInit[0] + P;

			// adjust temperature
			// liquid: must account for latent heat of fusion
			T[0] = (P *(T_air + LF/CI) + T[0] * mInit[0]) / mass;

			// adjust grid cell density
			d[0] = mass / dz[0];

			// if d > the density of ice, d = dIce
			if (d[0] > dIce){
				d[0] = dIce;           // adjust d
				dz[0] = mass / d[0];    // dz is adjusted to conserve mass
			}
		}

		// check for conservation of mass
		mass=0; for(int i=0;i<m;i++)mass+=d[i]*dz[i]; 

		mass_diff = mass - massinit - P;
		
		#ifndef _HAVE_ADOLC_  //avoid round operation. only check in forward mode.
		mass_diff = round(mass_diff * 100)/100;
		if (mass_diff > 0) _error_("mass not conserved in accumulation function");
		#endif

	}
	/*Free ressources:*/
	if(mInit)xDelete<IssmDouble>(mInit);

	/*Assign output pointers:*/
	*pT=T;
	*pdz=dz;
	*pd=d;
	*pW=W;
	*pa=a;
	*pre=re;
	*pgdn=gdn;
	*pgsp=gsp;
	*pm=m;
} /*}}}*/
void melt(IssmDouble* pM, IssmDouble* pR, IssmDouble* pmAdd, IssmDouble* pdz_add, IssmDouble** pT, IssmDouble** pd, IssmDouble** pdz, IssmDouble** pW, IssmDouble** pa, IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, int* pn, IssmDouble dzMin, IssmDouble zMax, IssmDouble zMin, IssmDouble zTop, int sid){ /*{{{*/

	//// MELT ROUTINE

	// Description:
	// computes the quantity of meltwater due to snow temperature in excess of
	// 0 deg C, determines pore water content and adjusts grid spacing

	/*intermediary:*/
	IssmDouble* m=NULL;
	IssmDouble* maxF=NULL;
	IssmDouble* dW=NULL;
	IssmDouble* exsW=NULL;
	IssmDouble* exsT=NULL;
	IssmDouble* surpT=NULL;
	IssmDouble* surpE=NULL;
	IssmDouble* F=NULL;
	IssmDouble* flxDn=NULL;
	IssmDouble  ER=0;
	IssmDouble* EI=NULL;
	IssmDouble* EW=NULL;
	IssmDouble* M=NULL;
	int*       D=NULL;
	
	IssmDouble sumM;
	IssmDouble sumER;
	IssmDouble addE;
	IssmDouble mSum0;
	IssmDouble sumE0;
	IssmDouble mSum1;
	IssmDouble sumE1;
	IssmDouble dE;
	IssmDouble dm;
	IssmDouble X;
	IssmDouble Wi;
    
	IssmDouble Ztot;
	IssmDouble T_bot;
	IssmDouble m_bot;
	IssmDouble dz_bot;
	IssmDouble d_bot;
	IssmDouble W_bot;
	IssmDouble a_bot;
	IssmDouble re_bot;
	IssmDouble gdn_bot;
	IssmDouble gsp_bot;
	bool        top=false;
    
	IssmDouble* Zcum=NULL;
	IssmDouble* dzMin2=NULL;
	IssmDouble zY2=1.025;
	bool lastCellFlag;
	int X1=0;
	int X2=0;
    
	int        D_size;

	/*outputs:*/
	IssmDouble  mAdd;
	IssmDouble dz_add;
	IssmDouble  Rsum;
	IssmDouble* T=*pT;
	IssmDouble* d=*pd;
	IssmDouble* dz=*pdz;
	IssmDouble* W=*pW;
	IssmDouble* a=*pa;
	IssmDouble* re=*pre;
	IssmDouble* gdn=*pgdn;
	IssmDouble* gsp=*pgsp;
	int         n=*pn;
	IssmDouble* R=0;
	
	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   melt module\n");

	//// INITIALIZATION

	/*Allocations: */
	M=xNewZeroInit<IssmDouble>(n); 
	maxF=xNew<IssmDouble>(n); 
	dW=xNew<IssmDouble>(n); 

	// specify constants
	const IssmDouble CtoK = 273.15;  // clecius to Kelvin conversion
	const IssmDouble CI = 2102;      // specific heat capacity of snow/ice (J kg-1 k-1)
	const IssmDouble LF = 0.3345E6;  // latent heat of fusion(J kg-1)
	const IssmDouble dPHC = 830;     // pore hole close off density[kg m-3]
	const IssmDouble dIce = 910;     // density of ice [kg m-3]

	// store initial mass [kg] and energy [J]
	m=xNew<IssmDouble>(n); for(int i=0;i<n;i++) m[i] = dz[i]* d[i];                    // grid cell mass [kg]
	EI=xNew<IssmDouble>(n); for(int i=0;i<n;i++)EI[i] = m[i] * T[i] * CI;               // initial enegy of snow/ice
	EW=xNew<IssmDouble>(n); for(int i=0;i<n;i++)EW[i]= W[i] * (LF + CtoK * CI);     // initial enegy of water

	mSum0 = cellsum(W,n) + cellsum(m,n);        // total mass [kg]
	sumE0 = cellsum(EI,n) + cellsum(EW,n);      // total energy [J]

	// initialize melt and runoff scalars
	Rsum = 0;       // runoff [kg]
	sumM = 0;       // total melt [kg]
	mAdd = 0;       // mass added/removed to/from base of model [kg]
	addE = 0;       // energy added/removed to/from base of model [J]
	dz_add=0;      // thickness of the layer added/removed to/from base of model [m]

	// calculate temperature excess above 0 deg C
	exsT=xNewZeroInit<IssmDouble>(n);
	for(int i=0;i<n;i++) exsT[i]= fmax(0, T[i] - CtoK);        // [K] to [°C]

	// new grid point center temperature, T [K]
	for(int i=0;i<n;i++) T[i]-=exsT[i];

	// specify irreducible water content saturation [fraction]
	const IssmDouble Swi = 0.07;                     // assumed constant after Colbeck, 1974

	//// REFREEZE PORE WATER
	// check if any pore water
	if (cellsum(W,n) > 0){
		if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("      pore water refreeze\n");
		// calculate maximum freeze amount, maxF [kg]
		for(int i=0;i<n;i++) maxF[i] = fmax(0, -((T[i] - CtoK) * m[i] * CI) / LF);

		// freeze pore water and change snow/ice properties
		for(int i=0;i<n;i++) dW[i] = fmin(maxF[i], W[i]);    // freeze mass [kg]   
		for(int i=0;i<n;i++) W[i] -= dW[i];                                            // pore water mass [kg]
		for(int i=0;i<n;i++) m[i] += dW[i];                                            // new mass [kg]
		for(int i=0;i<n;i++) d[i] = m[i] / dz[i];                                    // density [kg m-3]   
		for(int i=0;i<n;i++) T[i] = T[i] + (dW[i]*(LF+(CtoK - T[i])*CI)/(m[i]*CI));      // temperature [K]

		// if pore water froze in ice then adjust d and dz thickness
		for(int i=0;i<n;i++)if(d[i]>dIce)d[i]=dIce;
		for(int i=0;i<n;i++) dz[i]= m[i]/d[i];
	}

	// squeeze water from snow pack
	exsW=xNew<IssmDouble>(n); 
	for(int i=0;i<n;i++){
		Wi= (910 - d[i]) * Swi * (m[i] / d[i]);        // irreducible water content [kg]
		exsW[i] = fmax(0, W[i] - Wi);                  // water "squeezed" from snow [kg]
	}

	//// MELT, PERCOLATION AND REFREEZE

	// run melt algorithm if there is melt water or excess pore water
	if ((cellsum(exsT,n) > 0) || (cellsum(exsW,n) > 0)){
		// _printf_(""MELT OCCURS");
		// check to see if thermal energy exceeds energy to melt entire cell
		// if so redistribute temperature to lower cells (temperature surplus)
		// (maximum T of snow before entire grid cell melts is a constant
		// LF/CI = 159.1342)
		surpT=xNew<IssmDouble>(n); for(int i=0;i<n;i++)surpT[i] = fmax(0, exsT [i]- 159.1342);
        
		if (cellsum(surpT,n) > 0 ){
			// _printf_("T Surplus");
			// calculate surplus energy
			surpE=xNew<IssmDouble>(n); for(int i=0;i<n;i++)surpE[i] = surpT[i] * CI * m[i];
            
			int i = 0;
			while (cellsum(surpE,n) > 0){
				// use surplus energy to increase the temperature of lower cell
				T[i+1] = surpE[i]/m[i+1]/CI + T[i+1];
                
				exsT[i+1] = fmax(0, T[i+1] - CtoK) + exsT[i+1];
				T[i+1] = fmin(CtoK, T[i+1]);
                
				surpT[i+1] = fmax(0, exsT[i+1] - 159.1342);
				surpE[i+1] = surpT[i+1] * CI * m[i+1];
                
				// adjust current cell properties (again 159.1342 is the max T)
				exsT[i] = 159.1342;
				surpE[i] = 0;
				i = i + 1;
			}
		}

		// convert temperature excess to melt [kg]
		for(int i=0;i<n;i++) M[i] = exsT[i] * d[i] * dz[i] * CI / LF;      // melt
		sumM = cellsum(M,n);                                               // total melt [kg]

		// calculate maximum refreeze amount, maxF [kg]
		for(int i=0;i<n;i++)maxF[i] = fmax(0, -((T[i] - CtoK) * d[i] * dz[i] * CI)/ LF);

		// initialize refreeze, runoff, flxDn and dW vectors [kg]
 		IssmDouble* F = xNewZeroInit<IssmDouble>(n);
		IssmDouble* R=xNewZeroInit<IssmDouble>(n);

		for(int i=0;i<n;i++)dW[i] = 0;
		flxDn=xNewZeroInit<IssmDouble>(n+1); for(int i=0;i<n;i++)flxDn[i+1]=F[i];

		// determine the deepest grid cell where melt/pore water is generated
		X = 0;
		for(int i=n-1;i>=0;i--){
			if(M[i]>0 || reCast<int,IssmDouble>(exsW[i])){
				X=i;
				break;
			}
		}

		//// meltwater percolation
		for(int i=0;i<n;i++){
			// calculate total melt water entering cell
			IssmDouble inM = M[i]+ flxDn[i];

			// break loop if there is no meltwater and if depth is > mw_depth
			if (inM == 0 && i > X){
				break;
			}

			// if reaches impermeable ice layer all liquid water runs off (R)
			else if (d[i] >= dIce){  // dPHC = pore hole close off [kg m-3]
				// _printf_("ICE LAYER");
				// no water freezes in this cell
				// no water percolates to lower cell
				// cell ice temperature & density do not change

				m[i] = m[i] - M[i];                     // mass after melt
				Wi = (910-d[i]) * Swi * (m[i]/d[i]);    // irreducible water 
				dW[i] = fmin(inM, Wi - W[i]);            // change in pore water
				R[i] = fmax(0, inM - dW[i]);             // runoff
			}
			// check if no energy to refreeze meltwater
			else if (maxF[i] == 0){
				// _printf_("REFREEZE == 0");
				// no water freezes in this cell
				// cell ice temperature & density do not change

				m[i] = m[i] - M[i];                     // mass after melt
				Wi = (910-d[i]) * Swi * (m[i]/d[i]);    // irreducible water 
				dW[i] = fmin(inM, Wi-W[i]);              // change in pore water
				flxDn[i+1] = fmax(0, inM-dW[i]);         // meltwater out
				F[i] = 0;                               // no freeze 
			}
			// some or all meltwater refreezes
			else{
				// change in density density and temperature
				// _printf_("MELT REFREEZE");
				//-----------------------melt water-----------------------------
				IssmDouble dz_0 = m[i]/d[i];          
				IssmDouble dMax = (dIce - d[i])*dz_0;              // d max = dIce
				IssmDouble F1 = fmin(fmin(inM,dMax),maxF[i]);         // maximum refreeze               
				m[i] = m[i] + F1;                       // mass after refreeze
				d[i] = m[i]/dz_0;

				//-----------------------pore water-----------------------------
				Wi = (910-d[i])* Swi * dz_0;            // irreducible water 
				dW[i] = fmin(inM - F1, Wi-W[i]);         // change in pore water
				if (-dW[i]>W[i] ){
					dW[i]= W[i];
				}
				IssmDouble F2 = 0;                                 

				if (dW[i] < 0){                         // excess pore water
					dMax = (dIce - d[i])*dz_0;          // maximum refreeze                                             
					IssmDouble maxF2 = fmin(dMax, maxF[i]-F1);      // maximum refreeze
					F2 = fmin(-dW[i], maxF2);            // pore water refreeze
					m[i] = m[i] + F2;                   // mass after refreeze
					d[i] = m[i]/dz_0;
				}

				flxDn[i+1] = inM - F1 - dW[i] - F2;     // meltwater out        
				T[i] = T[i] + ((F1+F2)*(LF+(CtoK - T[i])*CI)/(m[i]*CI));// change in temperature


				// check if an ice layer forms 
				if (d[i] == dIce){
					// _printf_("ICE LAYER FORMS");
					// excess water runs off
					R[i] = flxDn[i+1];
					// no water percolates to lower cell
					flxDn[i+1] = 0;
				}
			}
		}


		//// GRID CELL SPACING AND MODEL DEPTH
		for(int i=0;i<n;i++)if (W[i] < 0) _error_("negative pore water generated in melt equations");
		
		// delete all cells with zero mass
		// adjust pore water
		for(int i=0;i<n;i++)W[i] += dW[i];

		// delete all cells with zero mass
		D_size=0; for(int i=0;i<n;i++)if(m[i]!=0)D_size++; D=xNew<int>(D_size);
		D_size=0; for(int i=0;i<n;i++)if(m[i]!=0){ D[D_size] = i; D_size++;}
		
		celldelete(&m,n,D,D_size);
		celldelete(&W,n,D,D_size);
		celldelete(&d,n,D,D_size);
		celldelete(&T,n,D,D_size);
		celldelete(&a,n,D,D_size);
		celldelete(&re,n,D,D_size);
		celldelete(&gdn,n,D,D_size);
		celldelete(&gsp,n,D,D_size);
		celldelete(&EI,n,D,D_size);
		celldelete(&EW,n,D,D_size);
		celldelete(&R,n,D,D_size);
		n=D_size;
		xDelete<int>(D);
	
		// calculate new grid lengths
		for(int i=0;i<n;i++)dz[i] = m[i] / d[i];

		//calculate Rsum:
		Rsum=cellsum(R,n);

		/*Free ressources:*/
		xDelete<IssmDouble>(F);
		xDelete<IssmDouble>(R);
	}
    
	//Merging of cells as they are burried under snow.
	Zcum=xNew<IssmDouble>(n);
	dzMin2=xNew<IssmDouble>(n);
    
	Zcum[0]=dz[0]; // Compute a cumulative depth vector
    
	for (int i=1;i<n;i++){
		Zcum[i]=Zcum[i-1]+dz[i];
	}
    
	for (int i=0;i<n;i++){
		if (Zcum[i]<=zTop){
			dzMin2[i]=dzMin;
		}
		else{
			dzMin2[i]=zY2*dzMin2[i-1];
		}
	}

	// check if depth is too small
	X = 0;
	for(int i=n-1;i>=0;i--){
		if(dz[i]<dzMin2[i]){
			X=i;
			break;
		}
	}
    
	//Last cell has to be treated separately if has to be merged (no underlying cell so need to merge with overlying cell)
	if(X==n-1){
		lastCellFlag = true;
		X = X-1;
	}
	else{
		lastCellFlag = false;
	}

	for (int i = 0; i<=X;i++){
		if (dz [i] < dzMin2[i]){                              // merge top cells
			// _printf_("dz > dzMin * 2')
			// adjust variables as a linearly weighted function of mass
			IssmDouble m_new = m[i] + m[i+1];
			T[i+1] = (T[i]*m[i] + T[i+1]*m[i+1]) / m_new;
			a[i+1] = (a[i]*m[i] + a[i+1]*m[i+1]) / m_new;
			re[i+1] = (re[i]*m[i] + re[i+1]*m[i+1]) / m_new;
			gdn[i+1] = (gdn[i]*m[i] + gdn[i+1]*m[i+1]) / m_new;
			gsp[i+1] = (gsp[i]*m[i] + gsp[i+1]*m[i+1]) / m_new;
            
			// merge with underlying grid cell and delete old cell
			dz [i+1] = dz[i] + dz[i+1];                 // combine cell depths
			d[i+1] = m_new / dz[i+1];                   // combine top densities
			W[i+1] = W[i+1] + W[i];                     // combine liquid water
			m[i+1] = m_new;                             // combine top masses
            
			// set cell to 99999 for deletion
			m[i] = 99999;
		}
	}

	//If last cell has to be merged
	if(lastCellFlag){
         //find closest cell to merge with
		for(int i=n-2;i>=0;i--){
			if(m[i]!=99999){
				X2=i;
				X1=n-1;
				break;
			}
		}
        
		// adjust variables as a linearly weighted function of mass
		IssmDouble m_new = m[X2] + m[X1];
		T[X1] = (T[X2]*m[X2] + T[X1]*m[X1]) / m_new;
		a[X1] = (a[X2]*m[X2] + a[X1]*m[X1]) / m_new;
		re[X1] = (re[X2]*m[X2] + re[X1]*m[X1]) / m_new;
		gdn[X1] = (gdn[X2]*m[X2] + gdn[X1]*m[X1]) / m_new;
		gsp[X1] = (gsp[X2]*m[X2] + gsp[X1]*m[X1]) / m_new;
        
		// merge with underlying grid cell and delete old cell
  		dz [X1] = dz[X2] + dz[X1];                 // combine cell depths
		d[X1] = m_new / dz[X1];                   // combine top densities
 		W[X1] = W[X1] + W[X2];                     // combine liquid water
		m[X1] = m_new;                             // combine top masses
        
		// set cell to 99999 for deletion
		m[X2] = 99999;
	}

	// delete combined cells
	D_size=0; for(int i=0;i<n;i++)if(m[i]!=99999)D_size++; D=xNew<int>(D_size); 
	D_size=0; for(int i=0;i<n;i++)if(m[i]!=99999){ D[D_size] = i; D_size++;}

	celldelete(&m,n,D,D_size);
	celldelete(&W,n,D,D_size);
	celldelete(&dz,n,D,D_size);
	celldelete(&d,n,D,D_size);
	celldelete(&T,n,D,D_size);
	celldelete(&a,n,D,D_size);
	celldelete(&re,n,D,D_size);
	celldelete(&gdn,n,D,D_size);
	celldelete(&gsp,n,D,D_size);
	celldelete(&EI,n,D,D_size);
	celldelete(&EW,n,D,D_size);
	n=D_size;
	xDelete<int>(D);
    
	// check if any of the top 10 cell depths are too large
	X=0;
	for(int i=9;i>=0;i--){
		if(dz[i]> 2* dzMin){
			X=i;
			break;
		}
	}
	
	int i=0;
	while(i<=X){
		if (dz [i] > dzMin *2){

				// _printf_("dz > dzMin * 2");
				// split in two
				cellsplit(&dz, n, i,.5);
				cellsplit(&W, n, i,.5);
				cellsplit(&m, n, i,.5);
				cellsplit(&T, n, i,1.0);
				cellsplit(&d, n, i,1.0);
				cellsplit(&a, n, i,1.0);
				cellsplit(&EI, n, i,1.0);
				cellsplit(&EW, n, i,1.0);
				cellsplit(&re, n, i,1.0);
				cellsplit(&gdn, n, i,1.0);
				cellsplit(&gsp, n, i,1.0);
				n++;
				X=X+1;
		}
		else i++;
	}

	//// CORRECT FOR TOTAL MODEL DEPTH
	// WORKS FINE BUT HAS BEEN DISABLED FOR CONVIENCE OF MODEL OUTPUT
	// INTERPRETATION
	// // calculate total model depth
	Ztot = cellsum(dz,n);
    
	if (Ztot < zMin){
		// printf("Total depth < zMin %f \n", Ztot);
		// mass and energy to be added
		mAdd = m[n-1]+W[n-1];
		addE = T[n-1]*m[n-1]*CI;
        
		// add a grid cell of the same size and temperature to the bottom
		dz_bot=dz[n-1];
		T_bot=T[n-1];
		W_bot=W[n-1];
		m_bot=m[n-1];
		d_bot=d[n-1];
		a_bot=a[n-1];
		re_bot=re[n-1];
		gdn_bot=gdn[n-1];
		gsp_bot=gsp[n-1];
        
		dz_add=dz_bot;
        
		newcell(&dz,dz_bot,top,n);
		newcell(&T,T_bot,top,n);
		newcell(&W,W_bot,top,n);
		newcell(&m,m_bot,top,n);
		newcell(&d,d_bot,top,n);
		newcell(&a,a_bot,top,n);
		newcell(&re,re_bot,top,n);
		newcell(&gdn,gdn_bot,top,n);
		newcell(&gsp,gsp_bot,top,n);
		n=n+1;
	}
	else if (Ztot > zMax){
		// printf("Total depth > zMax %f \n", Ztot);
		// mass and energy loss
		mAdd = -(m[n-1]+W[n-1]);
		addE = -(T[n-1]*m[n-1]*CI);
		dz_add=-(dz[n-1]);
        
		// add a grid cell of the same size and temperature to the bottom
		D_size=n-1;
		D=xNew<int>(D_size);
        
		for(int i=0;i<n-1;i++) D[i]=i;
		celldelete(&dz,n,D,D_size);
		celldelete(&T,n,D,D_size);
		celldelete(&W,n,D,D_size);
		celldelete(&m,n,D,D_size);
		celldelete(&d,n,D,D_size);
		celldelete(&a,n,D,D_size);
		celldelete(&re,n,D,D_size);
		celldelete(&gdn,n,D,D_size);
		celldelete(&gsp,n,D,D_size);
		n=D_size;
		xDelete<int>(D);
        
	}

	//// CHECK FOR MASS AND ENERGY CONSERVATION

	// calculate final mass [kg] and energy [J]
	sumER = Rsum * (LF + CtoK * CI);
	for(int i=0;i<n;i++)EI[i] = m[i] * T[i] * CI;
	for(int i=0;i<n;i++)EW[i] = W[i] * (LF + CtoK * CI);

	mSum1 = cellsum(W,n) + cellsum(m,n) + Rsum;
	sumE1 = cellsum(EI,n) + cellsum(EW,n);

	/*checks: */
	for(int i=0;i<n;i++) if (W[i]<0) _error_("negative pore water generated in melt equations\n");
	
	/*only in forward mode! avoid round in AD mode as it is not differentiable: */
	#ifndef _HAVE_ADOLC_
	dm = round(mSum0 - mSum1 + mAdd);
	dE = round(sumE0 - sumE1 - sumER +  addE);
	if (dm !=0  || dE !=0) _error_("mass or energy are not conserved in melt equations\n"
			<< "dm: " << dm << " dE: " << dE << "\n");
	#endif
	
	/*Free ressources:*/
	if(m)xDelete<IssmDouble>(m);
	if(EI)xDelete<IssmDouble>(EI);
	if(EW)xDelete<IssmDouble>(EW);
	if(maxF)xDelete<IssmDouble>(maxF);
	if(dW)xDelete<IssmDouble>(dW);
	if(exsW)xDelete<IssmDouble>(exsW);
	if(exsT)xDelete<IssmDouble>(exsT);
	if(surpT)xDelete<IssmDouble>(surpT);
	if(surpE)xDelete<IssmDouble>(surpE);
	if(flxDn)xDelete<IssmDouble>(flxDn);
	if(D)xDelete<int>(D);
	if(M)xDelete<IssmDouble>(M);
 	xDelete<IssmDouble>(Zcum);
	xDelete<IssmDouble>(dzMin2);
    
	/*Assign output pointers:*/
	*pM=sumM;
	*pR=Rsum;
	*pmAdd=mAdd;
	*pdz_add=dz_add;

	*pT=T;
	*pd=d;
	*pdz=dz;
	*pW=W;
	*pa=a;
	*pre=re;
	*pgdn=gdn;
	*pgsp=gsp;
	*pn=n;

} /*}}}*/ 
void densification(IssmDouble* d,IssmDouble* dz, IssmDouble* T, IssmDouble* re, int denIdx, IssmDouble C, IssmDouble dt, IssmDouble Tmean,IssmDouble dIce, int m, int sid){ /*{{{*/

	//// THIS NEEDS TO BE DOUBLE CHECKED AS THERE SEAMS TO BE LITTLE DENSIFICATION IN THE MODEL OUTOUT [MAYBE COMPATION IS COMPNSATED FOR BY TRACES OF SNOW???]

	//// FUNCTION INFO

	// Author: Alex Gardner, University of Alberta
	// Date last modified: FEB, 2008 

	// Description: 
	//   computes the densification of snow/firn using the emperical model of
	//   Herron and Langway (1980) or the semi-emperical model of Anthern et al.
	//   (2010)

	// Inputs:
	//   denIdx = densification model to use:
	//       1 = emperical model of Herron and Langway (1980)
	//       2 = semi-imerical model of Anthern et al. (2010)
	//       3 = physical model from Appendix B of Anthern et al. (2010)
	//   d   = initial snow/firn density [kg m-3]
	//   T   = temperature [K]
	//   dz  = grid cell size [m]
	//   C   = average accumulation rate [kg m-2 yr-1]
	//   dt  = time lapsed [s]
	//   re  = effective grain radius [mm];
	//   Ta  = mean annual temperature                                          

	// Reference: 
	// Herron and Langway (1980), Anthern et al. (2010)

	//// FOR TESTING
	// denIdx = 2;
	// d = 800;
	// T = 270;
	// dz = 0.005;
	// C = 200;
	// dt = 60*60;
	// re = 0.7;
	// Tmean = 273.15-18;

	//// MAIN FUNCTION
	// specify constants
	dt      = dt / 86400;  // convert from [s] to [d]
	// R     = 8.314        // gas constant [mol-1 K-1]
	// Ec    = 60           // activation energy for self-diffusion of water
	//                      // molecules through the ice tattice [kJ mol-1]
	// Eg    = 42.4         // activation energy for grain growth [kJ mol-1]

	/*intermediary: */
	IssmDouble c0,c1,H;
	
	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   densification module\n");

	// initial mass
	IssmDouble* mass_init = xNew<IssmDouble>(m);for(int i=0;i<m;i++) mass_init[i]=d[i] * dz[i];
	
	/*allocations and initialization of overburden pressure and factor H: */
	IssmDouble* cumdz = xNew<IssmDouble>(m-1);
	cumdz[0]=dz[0];
	for(int i=1;i<m-1;i++)cumdz[i]=cumdz[i-1]+dz[i];

	IssmDouble* obp = xNew<IssmDouble>(m);
	obp[0]=0;
	for(int i=1;i<m;i++)obp[i]=cumdz[i-1]*d[i-1];
	
	// calculate new snow/firn density for:
	//   snow with densities <= 550 [kg m-3]
	//   snow with densities > 550 [kg m-3]
		
	
	for(int i=0;i<m;i++){
		switch (denIdx){
			case 1: // Herron and Langway (1980)
				c0 = (11 * exp(-10160 / (T[i] * 8.314))) * C/1000;
				c1 = (575 * exp(-21400 / (T[i]* 8.314))) * pow(C/1000,.5);
				break;
			case 2: // Arthern et al. (2010) [semi-emperical]
				// common variable
				// NOTE: Ec=60000, Eg=42400 (i.e. should be in J not kJ)
				H = exp((-60000./(T[i] * 8.314)) + (42400./(Tmean * 8.314))) * (C * 9.81);
				c0 = 0.07 * H;
				c1 = 0.03 * H;
				break;

			case 3: // Arthern et al. (2010) [physical model eqn. B1]

				// common variable
				H = exp((-60/(T[i] * 8.314))) * obp[i] / pow(re[i]/1000,2.0);
				c0 = 9.2e-9 * H;
				c1 = 3.7e-9 * H;
				break;

			case 4: // Li and Zwally (2004)
				c0 = (C/dIce) * (139.21 - 0.542*Tmean)*8.36*pow(273.15 - T[i],-2.061);
				c1 = c0;
				break;

			case 5: // Helsen et al. (2008)
				// common variable
				c0 = (C/dIce) * (76.138 - 0.28965*Tmean)*8.36*pow(273.15 - T[i],-2.061);
				c1 = c0;
				break;
		}

		// new snow density
		if(d[i] <= 550) d[i] = d[i] + (c0 * (dIce - d[i]) / 365 * dt);
		else            d[i] = d[i] + (c1 * (dIce - d[i]) / 365 * dt);

		//disp((num2str(nanmean(c0 .* (dIce - d(idx)) / 365 * dt))))

		// do not allow densities to exceed the density of ice
		if(d[i]>dIce)d[i]=dIce;

		// calculate new grid cell length
		dz[i] = mass_init[i] / d[i];
	}

	/*Free ressources:*/
	xDelete<IssmDouble>(mass_init);
	xDelete<IssmDouble>(cumdz);
	xDelete<IssmDouble>(obp);

} /*}}}*/
void turbulentFlux(IssmDouble* pshf, IssmDouble* plhf, IssmDouble* pEC, IssmDouble Ta, IssmDouble Ts, IssmDouble V, IssmDouble eAir, IssmDouble pAir, IssmDouble ds, IssmDouble Ws, IssmDouble Vz, IssmDouble Tz, int sid){ /*{{{*/

	//// TURBULENT HEAT FLUX

	// Description: 
	// computed the surface sensible and latent heat fluxes [W m-2], this
	// function also calculated the mass loss/acreation due to
	// condensation/evaporation [kg]

	// Reference: 
	// Dingman, 2002.

	//// INPUTS:
	//   Ta: 2m air temperature [K]
	//   Ts: snow/firn/ice surface temperature [K]
	//   V: wind speed [m s^-^1]
	//   eAir: screen level vapor pressure [Pa]
	//   pAir: surface pressure [Pa]
	//   ds: surface density [kg/m^3]
	//   Ws: surface liquid water content [kg/m^2]
	//   Vz: height above ground at which wind (V) eas sampled [m]
	//   Tz: height above ground at which temperature (T) was sampled [m]

	//// FUNCTION INITILIZATION 

	// CA = 1005;                    // heat capacity of air (J kg-1 k-1)
	// LF = 0.3345E6;                // latent heat of fusion(J kg-1)
	// LV = 2.495E6;                 // latent heat of vaporization(J kg-1)
	// dIce = 910;                   // density of ice [kg m-3]
	
	/*intermediary:*/
	IssmDouble d_air;
	IssmDouble Ri;
	IssmDouble z0;
	IssmDouble coef_M,coef_H;
	IssmDouble An, C;
	IssmDouble L, eS;

	/*output: */
	IssmDouble shf, lhf, EC;
	
	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   turbulentFlux module\n");

	// calculated air density [kg/m3]
	d_air = 0.029 * pAir /(8.314 * Ta);

	//// Determine Surface Roughness
	// Bougamont, 2006
	// wind/temperature surface roughness height [m]
	if (ds < 910 && Ws == 0) z0 = 0.00012;               // 0.12 mm for dry snow
	else if (ds >= 910) z0 = 0.0032;                // 3.2 mm for ice 
	else z0 = 0.0013;                // 1.3 mm for wet snow

	//// Monin–Obukhov Stability Correction
	// Reference:
	// Ohmura, A., 1982: Climate and Energy-Balance on the Arctic Tundra.
	// Journal of Climatology, 2, 65-84.

	// if V = 0 goes to infinity therfore if V = 0 change
	if(V< .01) V=.01;

	// calculate the Bulk Richardson Number (Ri)
	Ri = (2*9.81* (Vz - z0) * (Ta - Ts)) / ((Ta + Ts)* pow(V,2));

	// calculate Monin–Obukhov stability factors 'coef_M' and 'coef_H'

	// do not allow Ri to exceed 0.19
	if(Ri>.19)Ri= 0.19;

	// calculate momentum 'coef_M' stability factor
	if (Ri > 0) coef_M = pow(1-5.2*Ri,-1); // if stable
	else coef_M = pow(1-18*Ri,-0.25);

	// calculate heat/wind 'coef_H' stability factor
	if (Ri < -0.03) coef_H = 1.3 * coef_M;
	else coef_H = coef_M;
		
	//// Bulk-transfer coefficient
	An =  pow(0.4,2) / pow(log(Tz/z0),2);     // Bulk-transfer coefficient
	C = An * d_air * V;             // shf & lhf common coefficient

	//// Sensible Heat
	// calculate the sensible heat flux [W m-2](Patterson, 1998)
	shf = C * 1005 * (Ta - Ts);

	// adjust using Monin–Obukhov stability theory
	shf = shf / (coef_M * coef_H);

	//// Latent Heat
	// determine if snow pack is melting & calcualte surface vapour pressure
	// over ice or liquid water
	if (Ts >= 273.15){
		L = 2.495E6;

		// for an ice surface Murphy and Koop, 2005 [Equation 7]
		eS = exp(9.550426 - 5723.265/Ts + 3.53068 * log(Ts)- 0.00728332 * Ts);
	}
	else{
		L = 2.8295E6; // latent heat of sublimation
		// for liquid surface (assume liquid on surface when Ts == 0 deg C)
		// Wright (1997), US Meteorological Handbook from Murphy and Koop,
		// 2005 Apendix A
		eS = 611.21 * exp(17.502 * (Ts - 273.15) / (240.97 + Ts - 273.15));
	}

	// Latent heat flux [W m-2]
	lhf = C * L * (eAir - eS) * 0.622 / pAir;

	// adjust using Monin–Obukhov stability theory (if lhf '+' then there is
	// energy and mass gained at the surface, if '-' then there is mass and 
	// energy loss at the surface. 
	lhf = lhf / (coef_M * coef_H);

	// mass loss (-)/acreation(+) due to evaporation/condensation [kg]
	EC = lhf * 86400 / L;

	/*assign output poitners: */
	*pshf=shf;
	*plhf=lhf;
	*pEC=EC;

} /*}}}*/
