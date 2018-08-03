import numpy as np
from TMeltingPoint  import TMeltingPoint

def DepthAvgTempCond(md):
   ''' compute conduction dependent temperature profile for an ice sheet. 
   Usage:
   Tbar=DepthAvgTempCond(md)
   '''

   Tpmp=TMeltingPoint(md.materials.meltingpoint,0) #pressure melting point at 0 pressure.

   k=md.materials.thermalconductivity
   G=md.basalforcings.geothermalflux
   H=md.geometry.thickness
   Ts=md.initialization.temperature
   alpha=G*H/k

   Tbar=np.zeros(md.mesh.numberofvertices,)

   #find temperature average when we are below melting point: 
   pos=np.nonzero( Ts+alpha < Tpmp)
   if pos:
	   Tbar[pos]=Ts[pos]+alpha[pos]/2 

   pos=np.nonzero( Ts+alpha>= Tpmp)
   if pos:
	   Tbar[pos]=Tpmp+(Tpmp**2-Ts[pos]**2)/2/alpha[pos]+ Tpmp*(Ts[pos]-Tpmp)/alpha[pos]
   
   #on ice shelf, easier: 
   pos=np.nonzero(md.mask.groundedice_levelset[0]<=0)
   if pos:
	   Tbar[pos]=(Ts[pos]+Tpmp)/2

   return Tbar
