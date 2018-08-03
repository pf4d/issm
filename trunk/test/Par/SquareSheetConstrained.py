import os.path
import numpy as np
import inspect
from verbose import verbose
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from paterson import paterson
from SetIceSheetBC import SetIceSheetBC
from arch import *

#Start defining model parameters here

#Geometry
hmin=300.
hmax=1000.
ymin=np.min(md.mesh.y)
ymax=np.max(md.mesh.y)
xmin=np.min(md.mesh.x)
xmax=np.max(md.mesh.x)
md.geometry.thickness=hmax+(hmin-hmax)*(md.mesh.y-ymin)/(ymax-ymin)+0.1*(hmin-hmax)*(md.mesh.x-xmin)/(xmax-xmin)
md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness+20.
md.geometry.surface=md.geometry.base+md.geometry.thickness

#Initial velocity 
x         = np.array(archread('../Data/SquareSheetConstrained.arch','x'))
y         = np.array(archread('../Data/SquareSheetConstrained.arch','y'))
vx        = np.array(archread('../Data/SquareSheetConstrained.arch','vx'));
vy        = np.array(archread('../Data/SquareSheetConstrained.arch','vy'));
index     = archread('../Data/SquareSheetConstrained.arch','index').astype(int);

md.initialization.vx=InterpFromMeshToMesh2d(index,x,y,vx,md.mesh.x,md.mesh.y)[0]
md.initialization.vy=InterpFromMeshToMesh2d(index,x,y,vy,md.mesh.x,md.mesh.y)[0]
md.initialization.vz=np.zeros((md.mesh.numberofvertices))
md.initialization.pressure=np.zeros((md.mesh.numberofvertices))

#Materials
md.initialization.temperature=(273.-20.)*np.ones((md.mesh.numberofvertices))
md.materials.rheology_B=paterson(md.initialization.temperature)
md.materials.rheology_n=3.*np.ones((md.mesh.numberofelements))

#Calving
md.calving.calvingrate=np.zeros((md.mesh.numberofvertices))
md.levelset.spclevelset=np.nan*np.ones((md.mesh.numberofvertices))

#Friction
md.friction.coefficient=20.*np.ones((md.mesh.numberofvertices))
md.friction.coefficient[np.where(md.mask.groundedice_levelset<0.)[0]]=0.
md.friction.p=np.ones((md.mesh.numberofelements))
md.friction.q=np.ones((md.mesh.numberofelements))

#Numerical parameters
md.stressbalance.viscosity_overshoot=0.0
md.masstransport.stabilization=1.
md.thermal.stabilization=1.
md.verbose=verbose(0)
md.settings.waitonlock=30
md.stressbalance.restol=0.05
md.steadystate.reltol=0.05
md.stressbalance.reltol=0.05
md.stressbalance.abstol=np.nan
md.timestepping.time_step=1.
md.timestepping.final_time=3.

#GIA:
md.gia.lithosphere_thickness=100.*np.ones((md.mesh.numberofvertices)); # in km
md.gia.mantle_viscosity=1.*10**21*np.ones((md.mesh.numberofvertices)); # in Pa.s
md.materials.lithosphere_shear_modulus=6.7*10**10;                          # in Pa
md.materials.lithosphere_density=3.32;                                      # in g/cm^-3
md.materials.mantle_shear_modulus=1.45*10**11;                              # in Pa
md.materials.mantle_density=3.34;                                           # in g/cm^-3

#Boundary conditions:
md=SetIceSheetBC(md)

#Change name so that no test have the same name
if len(inspect.stack()) > 2:
	md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
