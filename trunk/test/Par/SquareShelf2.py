import os.path
import inspect
from arch import *
import numpy
from verbose import verbose
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from paterson import paterson
from SetIceShelfBC import SetIceShelfBC

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
md.geometry.surface=md.geometry.base+md.geometry.thickness

#Initial velocity and pressure
x         = numpy.array(archread('../Data/SquareShelf.arch','x'))
y         = numpy.array(archread('../Data/SquareShelf.arch','y'))
vx        = numpy.array(archread('../Data/SquareShelf.arch','vx'));
vy        = numpy.array(archread('../Data/SquareShelf.arch','vy'));
index     = archread('../Data/SquareShelf.arch','index').astype(int);
#dbg - begin
# #print 'vars in SquareShelf.nc:'
# #for v in iVelF.variables:
# #	print v
#dbg - end 

[md.initialization.vx]=InterpFromMeshToMesh2d(index,x,y,vx,md.mesh.x,md.mesh.y)
[md.initialization.vy]=InterpFromMeshToMesh2d(index,x,y,vy,md.mesh.x,md.mesh.y)
md.initialization.vz=numpy.zeros((md.mesh.numberofvertices))
md.initialization.pressure=numpy.zeros((md.mesh.numberofvertices))

#dbg - begin
#print '...vx:'
#print md.initialization.vx
#print '...vy:'
#print md.initialization.vy
##print '...vz:'
##print md.initialization.vz
##print '...pressure:'
##print md.initialization.pressure
#dbg - end 


#Materials
md.initialization.temperature = (273.-20.)*numpy.ones((md.mesh.numberofvertices))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3.*numpy.ones((md.mesh.numberofelements))

#Friction
md.friction.coefficient = 20.*numpy.ones((md.mesh.numberofvertices))
md.friction.coefficient[numpy.nonzero(md.mask.groundedice_levelset<0.)[0]]=0.
md.friction.p = numpy.ones((md.mesh.numberofelements))
md.friction.q = numpy.ones((md.mesh.numberofelements))

#Numerical parameters
md.stressbalance.viscosity_overshoot = 0.3
md.masstransport.stabilization = 1.
md.thermal.stabilization = 1.
md.settings.waitonlock = 30
md.verbose=verbose()
md.stressbalance.restol = 0.10
md.steadystate.reltol = 0.02
md.stressbalance.reltol = 0.02
md.stressbalance.abstol = float('nan')
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.

#Boundary conditions:
# #md=SetIceShelfBC(md)
md=SetIceShelfBC(md,'../Exp/SquareFront2.exp')

#Change name so that no test have the same name
if len(inspect.stack()) > 2:
	md.miscellaneous.name=os.path.basename(inspect.stack()[2][1]).split('.')[0]
