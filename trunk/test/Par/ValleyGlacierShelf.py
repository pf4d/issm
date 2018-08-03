import os.path
from arch import *
import numpy
import inspect
import math
from verbose import verbose
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from paterson import paterson
from SetIceShelfBC import SetIceShelfBC

#Start defining model parameters here
x=md.mesh.x
y=md.mesh.y
xmin, xmax = min(x), max(x)
ymin, ymax = min(y), max(y)
Lx=(xmax-xmin)
Ly=(ymax-ymin)
xm,ym = (xmin+xmax)/2., (ymin+ymax)/2.

#Geometry: U-shaped valley in y direction
thk_center = 1000.
thk_margin = 0.5*thk_center
bmax=0.
bmin=-thk_center*md.materials.rho_ice/md.materials.rho_water

alpha=2./3.
slope = 0.9*(bmin-bmax)*(x-xmin)/(Lx*alpha) + 0.1*(bmin-bmax)*(y-ymin)/(Ly) + bmax
md.geometry.surface= (thk_center+bmax) + slope 
md.geometry.base=bmax + slope + 4./Ly**2*(thk_center-thk_margin)*(numpy.power(y-ym,2))
md.geometry.thickness=md.geometry.surface - md.geometry.base
md.geometry.bed = md.geometry.base

#Mask
md.mask.ice_levelset=x - alpha*Lx
md.mask.groundedice_levelset= numpy.ones((md.mesh.numberofvertices))

#Initial velocity 
md.initialization.vx=numpy.zeros((md.mesh.numberofvertices))
md.initialization.vy=numpy.zeros((md.mesh.numberofvertices))
md.initialization.vz=numpy.zeros((md.mesh.numberofvertices))
md.initialization.pressure=numpy.zeros((md.mesh.numberofvertices))

#Materials
md.initialization.temperature=(273.15-5.)*numpy.ones((md.mesh.numberofvertices))
md.initialization.waterfraction=numpy.zeros((md.mesh.numberofvertices))
md.initialization.watercolumn=numpy.zeros((md.mesh.numberofvertices))
md.materials.rheology_B=paterson(md.initialization.temperature)
md.materials.rheology_n=3.*numpy.ones((md.mesh.numberofelements))

#Thermal
md.thermal.isenthalpy=False
md.thermal.spctemperature=float('nan')*numpy.ones((md.mesh.numberofvertices))

#Groundingline
md.groundingline.migration='SubelementMigration'

#Surface mass balance and basal melting
md.smb.mass_balance=0.3*numpy.ones((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate=md.smb.mass_balance
md.basalforcings.floatingice_melting_rate=md.smb.mass_balance

#Friction
md.friction.coefficient=20.*numpy.ones((md.mesh.numberofvertices))
md.friction.coefficient[numpy.nonzero(md.mask.groundedice_levelset<0.)[0]]=0.
md.friction.p=numpy.ones((md.mesh.numberofelements))
md.friction.q=numpy.ones((md.mesh.numberofelements))

#Transient
md.transient.isstressbalance=True
md.transient.ismovingfront=True
md.transient.ismasstransport=False
md.transient.isthermal=False
md.transient.isgroundingline=True
md.transient.isgia=False

#Stressbalance
md.stressbalance.maxiter=100
md.stressbalance.viscosity_overshoot=0.0
md.stressbalance.restol=0.05
md.stressbalance.reltol=0.05
md.stressbalance.abstol=float('nan')

#Masstransport
md.calving.calvingrate=0.*numpy.ones((md.mesh.numberofvertices))
md.calving.meltingrate=0.*numpy.ones((md.mesh.numberofvertices))
md.levelset.spclevelset=float('NaN')*numpy.ones((md.mesh.numberofvertices))
md.masstransport.stabilization=1.

#Numerical parameters
md.thermal.stabilization=1.
md.settings.waitonlock=30
md.steadystate.reltol=0.05
md.timestepping.time_step=1.
md.timestepping.final_time=3.

#Verbose
md.verbose = verbose(0)

#Deal with boundary conditions:
md = SetIceShelfBC(md)

#Change name so that no tests have the same name
if len(inspect.stack()) > 2:
	md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
