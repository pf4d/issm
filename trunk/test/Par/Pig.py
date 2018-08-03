import os.path
import inspect
from arch import *
import numpy
from verbose import verbose
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from paterson import paterson 
from SetMarineIceSheetBC import SetMarineIceSheetBC

#Start defining model parameters here

#Geometry and observation
x         = numpy.array(archread('../Data/Pig.arch','x'))
y         = numpy.array(archread('../Data/Pig.arch','y'))
vx_obs    = numpy.array(archread('../Data/Pig.arch','vx_obs')) 
vy_obs    = numpy.array(archread('../Data/Pig.arch','vy_obs')) 
index     = numpy.array(archread('../Data/Pig.arch','index')).astype(int)
surface   = numpy.array(archread('../Data/Pig.arch','surface'))
thickness = numpy.array(archread('../Data/Pig.arch','thickness')) 

[md.inversion.vx_obs]   =InterpFromMeshToMesh2d(index,x,y,vx_obs,md.mesh.x,md.mesh.y)
[md.inversion.vy_obs]   =InterpFromMeshToMesh2d(index,x,y,vy_obs,md.mesh.x,md.mesh.y)
[md.geometry.surface]  =InterpFromMeshToMesh2d(index,x,y,surface,md.mesh.x,md.mesh.y)
[md.geometry.thickness]=InterpFromMeshToMesh2d(index,x,y,thickness,md.mesh.x,md.mesh.y)
md.geometry.base=md.geometry.surface-md.geometry.thickness
md.initialization.vx=md.inversion.vx_obs
md.initialization.vy=md.inversion.vy_obs
md.initialization.vz=numpy.zeros((md.mesh.numberofvertices))
md.initialization.pressure=numpy.zeros((md.mesh.numberofvertices))

#Materials
md.initialization.temperature=(273.-20.)*numpy.ones((md.mesh.numberofvertices))
md.materials.rheology_B=paterson(md.initialization.temperature)
md.materials.rheology_n=3.*numpy.ones((md.mesh.numberofelements))
md.initialization.temperature=md.initialization.temperature

#Friction
md.friction.coefficient=50.*numpy.ones((md.mesh.numberofvertices))
md.friction.coefficient[numpy.nonzero(md.mask.groundedice_levelset<0.)[0]]=0.
md.friction.p=numpy.ones((md.mesh.numberofelements))
md.friction.q=numpy.ones((md.mesh.numberofelements))

#Numerical parameters
md.stressbalance.viscosity_overshoot=0.3
md.masstransport.stabilization=1.
md.verbose=verbose(0)
md.settings.waitonlock=30
md.timestepping.time_step=1.
md.timestepping.final_time=2.
md.stressbalance.restol=0.05
md.stressbalance.reltol=1.
md.steadystate.reltol=1.
md.stressbalance.abstol=float('nan')

#Boundary conditions:
md=SetMarineIceSheetBC(md)

#Change name so that no test have the same name
if len(inspect.stack()) > 2:
	md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
