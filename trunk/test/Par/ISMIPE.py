import numpy
from arch import *
from SetIceSheetBC import SetIceSheetBC

#Ok, start defining model parameters here

print "      creating thickness"
data = numpy.array(archread('../Data/ISMIPE.arch','data'));
md.geometry.surface=numpy.zeros((md.mesh.numberofvertices))
md.geometry.base=numpy.zeros((md.mesh.numberofvertices))
for i in xrange(0,md.mesh.numberofvertices):
	y=md.mesh.y[i]
	point1=numpy.floor(y/100.)
	point2=numpy.minimum(point1+1,50)
	coeff=(y-(point1-1.)*100.)/100.
	md.geometry.base[i]=(1.-coeff)*data[point1,1]+coeff*data[point2,1]
	md.geometry.surface[i]=(1.-coeff)*data[point1,2]+coeff*data[point2,2]
md.geometry.thickness=md.geometry.surface-md.geometry.base
md.geometry.thickness[numpy.nonzero(numpy.logical_not(md.geometry.thickness))]=0.01
md.geometry.base=md.geometry.surface-md.geometry.thickness

print "      creating drag"
md.friction.coefficient=numpy.zeros((md.mesh.numberofvertices))
md.friction.p=numpy.ones((md.mesh.numberofelements))
md.friction.q=numpy.ones((md.mesh.numberofelements))

print "      creating flow law parameter"
md.materials.rheology_B=6.8067*10**7*numpy.ones((md.mesh.numberofvertices))
md.materials.rheology_n=3.*numpy.ones((md.mesh.numberofelements))

print "      boundary conditions for stressbalance model:"
#Create node on boundary first (because we can not use mesh)
md=SetIceSheetBC(md)
