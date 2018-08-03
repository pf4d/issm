import numpy
import os
from collections import OrderedDict
from expwrite import expwrite
from triangle import triangle

def roundmesh(md,radius,resolution):
	"""
	ROUNDMESH - create an unstructured round mesh 

	   This script will generate a structured round mesh
	   - radius     : specifies the radius of the circle in meters
	   - resolution : specifies the resolution in meters

	   Usage:
	      md=roundmesh(md,radius,resolution)
	"""

	#First we have to create the domain outline 

	#Get number of points on the circle
	pointsonedge=numpy.floor((2.*numpy.pi*radius) / resolution)

	#Calculate the cartesians coordinates of the points
	x_list=numpy.ones(pointsonedge)
	y_list=numpy.ones(pointsonedge)
	theta=numpy.linspace(0.,2.*numpy.pi,num=pointsonedge,endpoint=False)
	x_list=roundsigfig(radius*x_list*numpy.cos(theta),12)
	y_list=roundsigfig(radius*y_list*numpy.sin(theta),12)
	A=OrderedDict()
	A['x']=[x_list]
	A['y']=[y_list]
	A['density']=1.
	expwrite(A,'RoundDomainOutline.exp')

	#Call Bamg
	md=triangle(md,'RoundDomainOutline.exp',resolution)
	#md=bamg(md,'domain','RoundDomainOutline.exp','hmin',resolution)

	#move the closest node to the center
	pos=numpy.argmin(md.mesh.x**2+md.mesh.y**2)
	md.mesh.x[pos]=0.
	md.mesh.y[pos]=0.

	#delete domain
	os.remove('RoundDomainOutline.exp')

	return md

def roundsigfig(x,n):

	digits=numpy.ceil(numpy.log10(numpy.abs(x)))
	x=x/10.**digits
	x=numpy.round(x,decimals=n)
	x=x*10.**digits

	pos=numpy.nonzero(numpy.isnan(x))
	x[pos]=0.

	return x

