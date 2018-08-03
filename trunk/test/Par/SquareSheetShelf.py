import os.path
import inspect
from arch import *
import numpy
from verbose import verbose
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from paterson import paterson 
from SetMarineIceSheetBC import SetMarineIceSheetBC

#Start defining model parameters here

#Geometry
hmin=300.
hmax=1000.
ymin=min(md.mesh.y)
ymax=max(md.mesh.y)
xmin=min(md.mesh.x)
xmax=max(md.mesh.x)
md.geometry.thickness=hmax+(hmin-hmax)*(md.mesh.y-ymin)/(ymax-ymin)+0.1*(hmin-hmax)*(md.mesh.x-xmin)/(xmax-xmin)
md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness
bed_sheet=-md.materials.rho_ice/md.materials.rho_water*(hmax+(hmin-hmax)*(ymax/2-ymin)/(ymax-ymin))
pos=numpy.nonzero(md.mesh.y<=ymax/2.)
md.geometry.base[pos]=bed_sheet
md.geometry.surface=md.geometry.base+md.geometry.thickness

#Initial velocity 
x         = numpy.array(archread('../Data/SquareSheetShelf.arch','x'))
y         = numpy.array(archread('../Data/SquareSheetShelf.arch','y'))
vx        = numpy.array(archread('../Data/SquareSheetShelf.arch','vx'));
vy        = numpy.array(archread('../Data/SquareSheetShelf.arch','vy'));
index     = numpy.array(archread('../Data/SquareSheetShelf.arch','index')).astype(int);

[md.initialization.vx]  = InterpFromMeshToMesh2d(index,x,y,vx,md.mesh.x,md.mesh.y)
[md.initialization.vy]  = InterpFromMeshToMesh2d(index,x,y,vy,md.mesh.x,md.mesh.y)
md.initialization.vz=numpy.zeros((md.mesh.numberofvertices))
md.initialization.pressure=numpy.zeros((md.mesh.numberofvertices))

#Materials
md.initialization.temperature=(273.-20.)*numpy.ones((md.mesh.numberofvertices))
md.materials.rheology_B=paterson(md.initialization.temperature)
md.materials.rheology_n=3.*numpy.ones((md.mesh.numberofelements))

#Accumulation and melting
md.smb.mass_balance=10.*numpy.ones((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate=5.*numpy.ones((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate=5.*numpy.ones((md.mesh.numberofvertices))

#Friction
md.friction.coefficient=20.*numpy.ones((md.mesh.numberofvertices))
md.friction.coefficient[numpy.nonzero(md.mask.groundedice_levelset<0.)[0]]=0.
md.friction.p=numpy.ones((md.mesh.numberofelements))
md.friction.q=numpy.ones((md.mesh.numberofelements))

#Numerical parameters
md.stressbalance.viscosity_overshoot=0.0
md.masstransport.stabilization=1
md.thermal.stabilization=1
md.verbose=verbose(0)
md.settings.waitonlock=30
md.stressbalance.restol=0.05
md.steadystate.reltol=0.05
md.stressbalance.reltol=0.05
md.stressbalance.abstol=float('NaN')
md.timestepping.time_step=1.
md.timestepping.final_time=3.

#Deal with boundary conditions:
md=SetMarineIceSheetBC(md,'../Exp/SquareFront.exp')

#Change name so that no test have the same name
if len(inspect.stack()) > 2:
	md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
