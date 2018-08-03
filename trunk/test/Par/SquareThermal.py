import numpy
from paterson import paterson
from SetMarineIceSheetBC import SetMarineIceSheetBC

#Ok, start defining model parameters here

md.timestepping.time_step=0

print "      creating thickness"
h=1000.
md.geometry.thickness=h*numpy.ones((md.mesh.numberofvertices))
md.geometry.base=-1000.*numpy.ones((md.mesh.numberofvertices))
md.geometry.surface=md.geometry.base+md.geometry.thickness;

print "      creating velocities"
md.initialization.vx=numpy.zeros((md.mesh.numberofvertices))
md.initialization.vy=numpy.zeros((md.mesh.numberofvertices))
md.initialization.vz=numpy.zeros((md.mesh.numberofvertices))

print "      creating drag"
md.friction.coefficient=200.*numpy.ones((md.mesh.numberofvertices))
md.friction.coefficient[numpy.nonzero(md.mask.groundedice_levelset<0.)[0]]=0.
md.friction.p=numpy.ones((md.mesh.numberofelements))
md.friction.q=numpy.ones((md.mesh.numberofelements))

print "      creating temperatures"
md.initialization.temperature=(273.-20.)*numpy.ones((md.mesh.numberofvertices))

print "      creating flow law parameter"
md.materials.rheology_B=paterson(md.initialization.temperature)
md.materials.rheology_n=3.*numpy.ones((md.mesh.numberofelements))

print "      creating surface mass balance"
md.smb.mass_balance=numpy.ones((md.mesh.numberofvertices))/md.constants.yts    #1m/a
md.basalforcings.melting_rate=0.*numpy.ones((md.mesh.numberofvertices))/md.constants.yts    #1m/a

#Deal with boundary conditions:

print "      boundary conditions for stressbalance model"
md=SetMarineIceSheetBC(md,'../Exp/SquareFront.exp')

print "      boundary conditions for thermal model"
md.thermal.spctemperature=md.initialization.temperature
md.basalforcings.geothermalflux=numpy.zeros((md.mesh.numberofvertices)) 
md.basalforcings.geothermalflux[numpy.nonzero(md.mask.groundedice_levelset>0.)[0]]=1.*10**-3    #1 mW/m^2
