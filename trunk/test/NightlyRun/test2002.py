#Test Name: EarthSlr
from model import *
from socket import gethostname
import numpy as np
from parameterize import *
from solve import *
from gmshplanet import *
from maskpsl import *
from gmtmask import *
from paterson import *
from love_numbers import *

#mesh earth:
md=model()
md.mesh=gmshplanet('radius',6.371012*10**3,'resolution',700.) #500 km resolution mesh

#parameterize slr solution:
#slr loading:
md.slr.deltathickness=np.zeros((md.mesh.numberofelements))
md.slr.sealevel=np.zeros((md.mesh.numberofvertices))
#antarctica
late=np.sum(md.mesh.lat[md.mesh.elements-1],axis=1)/3
longe=np.sum(md.mesh.long[md.mesh.elements-1],axis=1)/3
pos=np.where(late <-80)
md.slr.deltathickness[pos]=-100
#greenland 
pos=np.where(np.logical_and.reduce((late > 70,late < 80,longe>-60,longe<-30)))
md.slr.deltathickness[pos]=-100

#elastic loading from love numbers:
nlov=101
md.slr.love_h = love_numbers('h')[:nlov]
md.slr.love_k = love_numbers('k')[:nlov]
md.slr.love_l = love_numbers('l')[:nlov]

#mask:
md.mask=maskpsl() # use maskpsl class (instead of mask) to store the ocean function as a ocean_levelset
mask=gmtmask(md.mesh.lat,md.mesh.long)

icemask=np.ones((md.mesh.numberofvertices))
pos=np.where(mask==0)[0]  
icemask[pos]=-1
pos=np.where(np.sum(mask[md.mesh.elements.astype(int)-1],axis=1)<3)[0]
icemask[md.mesh.elements[pos,:].astype(int)-1]=-1
md.mask.ice_levelset=icemask

md.mask.ocean_levelset=np.zeros((md.mesh.numberofvertices))
pos=np.where(md.mask.ice_levelset==1)
md.mask.ocean_levelset[pos]=1

#make sure that the ice level set is all inclusive:
md.mask.land_levelset=np.zeros((md.mesh.numberofvertices))
md.mask.groundedice_levelset=-np.ones((md.mesh.numberofvertices))

#make sure wherever there is an ice load, that the mask is set to ice:
pos=np.nonzero(md.slr.deltathickness)[0]
icemask[md.mesh.elements[pos,:]-1]=-1


#geometry
di=md.materials.rho_ice/md.materials.rho_water
md.geometry.thickness=np.ones((md.mesh.numberofvertices))
md.geometry.surface=(1-di)*np.zeros((md.mesh.numberofvertices))
md.geometry.base=md.geometry.surface-md.geometry.thickness
md.geometry.bed=md.geometry.base

#materials
md.initialization.temperature=273.25*np.ones((md.mesh.numberofvertices))
md.materials.rheology_B=paterson(md.initialization.temperature)
md.materials.rheology_n=3*np.ones((md.mesh.numberofelements))

#Miscellaneous
md.miscellaneous.name='test2002'

#Solution parameters
md.slr.reltol=np.nan
md.slr.abstol=1e-3

#eustatic run: 
md.slr.rigid=0
md.slr.elastic=0
md=solve(md,'Sealevelrise')
Seustatic=md.results.SealevelriseSolution.Sealevel

#eustatic + rigid run: 
md.slr.rigid=1
md.slr.elastic=0
md=solve(md,'Sealevelrise')
Srigid=md.results.SealevelriseSolution.Sealevel

#eustatic + rigid + elastic run: 
md.slr.rigid=1
md.slr.elastic=1
md=solve(md,'Sealevelrise')
Selastic=md.results.SealevelriseSolution.Sealevel

#eustatic + rigid + elastic + rotation run: 
md.slr.rigid=1
md.slr.elastic=1
md.slr.rotation=1
md=solve(md,'Sealevelrise')
Srotation=md.results.SealevelriseSolution.Sealevel

#Fields and tolerances to track changes
field_names     =['Eustatic','Rigid','Elastic','Rotation']
field_tolerances=[1e-13,1e-13,1e-13,1e-13]
field_values=[Seustatic,Srigid,Selastic,Srotation]
