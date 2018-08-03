import numpy
from SetIceSheetBC import SetIceSheetBC

#Ok, start defining model parameters here
md.verbose=2

print "      creating thickness"
md.geometry.surface=-md.mesh.x*numpy.tan(3.*numpy.pi/180.)
#md.geometry.base=md.geometry.surface-1000.
md.geometry.base=md.geometry.surface-1000.+100.*numpy.exp(-((md.mesh.x-numpy.max(md.mesh.x)/2.)**2+(md.mesh.y-numpy.max(md.mesh.y)/2.)**2)/(10000.**2))
md.geometry.thickness=md.geometry.surface-md.geometry.base

print "      creating drag"
md.friction.coefficient=numpy.sqrt(md.constants.yts/(2.140373*10**-7*1000.))*numpy.ones((md.mesh.numberofvertices))
md.friction.p=numpy.ones((md.mesh.numberofelements))
md.friction.q=numpy.zeros((md.mesh.numberofelements))

print "      creating flow law parameter"
md.materials.rheology_B=1.4734*10**14*numpy.ones((md.mesh.numberofvertices))
md.materials.rheology_n=1.*numpy.ones((md.mesh.numberofelements))
md.materials.rheology_law='None'

print "      boundary conditions for stressbalance model"
#Create node on boundary first (because we cannot use mesh)
md=SetIceSheetBC(md)
md.stressbalance.spcvx=100.*numpy.ones((md.mesh.numberofvertices))
md.initialization.vx=numpy.zeros((md.mesh.numberofvertices))
md.initialization.vy=numpy.zeros((md.mesh.numberofvertices))
md.initialization.vz=numpy.zeros((md.mesh.numberofvertices))
md.initialization.vel=numpy.zeros((md.mesh.numberofvertices))
md.initialization.pressure=numpy.zeros((md.mesh.numberofvertices))
md.initialization.temperature=255.*numpy.ones((md.mesh.numberofvertices))
pos=numpy.nonzero(numpy.logical_or(numpy.logical_or(md.mesh.x==numpy.min(md.mesh.x),md.mesh.x==numpy.max(md.mesh.x)),numpy.logical_or(md.mesh.y==numpy.min(md.mesh.y),md.mesh.y==numpy.max(md.mesh.y))))
md.balancethickness.spcthickness=float('NaN')*numpy.ones((md.mesh.numberofvertices))
md.balancethickness.spcthickness[pos]=md.geometry.thickness[pos]
md.masstransport.spcthickness=float('NaN')*numpy.ones((md.mesh.numberofvertices))
md.masstransport.spcthickness[pos]=md.geometry.thickness[pos]
md.thermal.spctemperature=255.*numpy.ones((md.mesh.numberofvertices))
md.basalforcings.geothermalflux=0.4*numpy.ones((md.mesh.numberofvertices))

#Parallel options
md.mesh.average_vertex_connectivity=200

#Transient options
md.timestepping.time_step=1.
md.timestepping.final_time=10.
md.masstransport.stabilization=1
md.thermal.stabilization=1
md.thermal.penalty_threshold=10**5
md.transient.isthermal=0
