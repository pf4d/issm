import numpy
from SetIceSheetBC import SetIceSheetBC

#Ok, start defining model parameters here

print "      creating thickness"
md.geometry.surface=-md.mesh.x*numpy.tan(0.5*numpy.pi/180.)
md.geometry.base=md.geometry.surface-1000.+500.*numpy.sin(md.mesh.x*2.*numpy.pi/numpy.max(md.mesh.x))
md.geometry.thickness=md.geometry.surface-md.geometry.base

print "      creating drag"
md.friction.coefficient=200.*numpy.ones((md.mesh.numberofvertices))
md.friction.coefficient[numpy.nonzero(md.mask.groundedice_levelset<0.)[0]]=0.
md.friction.p=numpy.ones((md.mesh.numberofelements))
md.friction.q=numpy.ones((md.mesh.numberofelements))

print "      creating flow law parameter"
md.materials.rheology_B=6.8067*10**7*numpy.ones((md.mesh.numberofvertices))
md.materials.rheology_n=3.*numpy.ones((md.mesh.numberofelements))

print "      boundary conditions for stressbalance model"
#Create node on boundary first (because we cannot use mesh)
md=SetIceSheetBC(md)
