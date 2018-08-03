/*!\file:  SurfaceMassBalancex.h
 * \brief header file for SMB
 */ 

#ifndef _SurfaceMassBalancex_H
#define _SurfaceMassBalancex_H

#include "../../classes/classes.h"

/* local prototypes: */
void SurfaceMassBalancex(FemModel* femmodel);
void SmbGradientsx(FemModel* femmodel);
void SmbGradientsElax(FemModel* femmodel);
void Delta18oParameterizationx(FemModel* femmodel);
void MungsmtpParameterizationx(FemModel* femmodel);
void Delta18opdParameterizationx(FemModel* femmodel);
void PositiveDegreeDayx(FemModel* femmodel);
void SmbHenningx(FemModel* femmodel);
void SmbComponentsx(FemModel* femmodel);
void SmbMeltComponentsx(FemModel* femmodel); 

/*GEMB: */
void       Gembx(FemModel* femmodel);
void       GembgridInitialize(IssmDouble** pdz, int* psize, IssmDouble zTop, IssmDouble dzTop, IssmDouble zMax, IssmDouble zY); 
IssmDouble Marbouty(IssmDouble T, IssmDouble d, IssmDouble dT);
void grainGrowth(IssmDouble* pre, IssmDouble* pgdn, IssmDouble* pgsp, IssmDouble* T,IssmDouble* dz,IssmDouble* d, IssmDouble* W,IssmDouble smb_dt,int m,int aIdx, int sid);
void albedo(IssmDouble* a,int aIdx, IssmDouble* re, IssmDouble* d, IssmDouble cldFrac, IssmDouble aIce, IssmDouble aSnow, IssmDouble* T, IssmDouble* W, IssmDouble P, IssmDouble EC, IssmDouble t0wet, IssmDouble t0dry, IssmDouble K, IssmDouble dt,int m, int sid);
void shortwave(IssmDouble** pswf, int swIdx, int aIdx, IssmDouble dsw, IssmDouble as, IssmDouble* d, IssmDouble* dz, IssmDouble* re, int m, int sid);
void thermo(IssmDouble* pEC, IssmDouble* T, IssmDouble* dz, IssmDouble* d, IssmDouble* swf, IssmDouble dlw, IssmDouble Ta, IssmDouble V, IssmDouble eAir, IssmDouble pAir, IssmDouble Ws, IssmDouble dt0, int m, IssmDouble Vz, IssmDouble Tz, int sid);
void accumulation(IssmDouble** pT, IssmDouble** pdz, IssmDouble** pd, IssmDouble** pW, IssmDouble** pa, IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, int* pm, IssmDouble Ta, IssmDouble P, IssmDouble dzMin, IssmDouble aSnow, int sid); 
void melt(IssmDouble* pM, IssmDouble* pR, IssmDouble* pmAdd, IssmDouble* pdz_add, IssmDouble** pT, IssmDouble** pd, IssmDouble** pdz, IssmDouble** pW, IssmDouble** pa, IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, int* pn, IssmDouble dzMin, IssmDouble zMax, IssmDouble zMin, IssmDouble zTop, int sid);
void densification(IssmDouble* d,IssmDouble* dz, IssmDouble* T, IssmDouble* re, int denIdx, IssmDouble C, IssmDouble dt, IssmDouble Tmean,IssmDouble rho_ice,int m, int sid);
void turbulentFlux(IssmDouble* pshf, IssmDouble* plhf, IssmDouble* pEC, IssmDouble Ta, IssmDouble Ts, IssmDouble V, IssmDouble eAir, IssmDouble pAir, IssmDouble ds, IssmDouble Ws, IssmDouble Vz, IssmDouble Tz, int sid);
#endif  /* _SurfaceMassBalancex_H*/
