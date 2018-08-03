from math import isnan
import MatlabFuncs as m
import numpy as np

def processmesh(md,data,options):
	"""
	PROCESSMESH - process the mesh for plotting

	Usage:
		x,y,z,elements,is2d=processmech(md,data,options)

	See also: PLOTMODEL, PROCESSDATA
	"""

	#some checks
	if md.mesh.numberofvertices==0:
		raise ValueError('processmesh error: mesh is empty')
	if md.mesh.numberofvertices==md.mesh.numberofelements:
		raise ValueError('processmesh error: the number of elements is the same as the number of nodes')

	if len(data)==0 or not isinstance(data,dict):
		
		if 'latlon' not in options.getfieldvalue('coord','xy').lower(): #convert to lower case for comparison
			x=md.mesh.x
			if 'x2d' in dir(md.mesh): x2d=md.mesh.x2d
			y=md.mesh.y
			if 'y2d' in dir(md.mesh): y2d=md.mesh.x2d
		else:
			x=md.mesh.long
			y=md.mesh.lat

		if 'z' in dir(md.mesh):
			z=md.mesh.z
		else:
			z=np.zeros_like(md.mesh.x)
		
		if 'elements2d' in dir(md.mesh): 
			elements2d=md.mesh.elements2d
			elements2d=elements2d-1  # subtract one since python indexes from zero
		elements=md.mesh.elements
		elements=elements-1

		#is it a 2D plot?
		if md.mesh.dimension()==2:
			is2d=1
		else:
			if options.getfieldvalue('layer',0)>=1:
				is2d=1
			else:
				is2d=0

		#layer projection?
		if options.getfieldvalue('layer',0)>=1:
			 if 'latlon' in options.getfieldvalue('coord','xy').lower():
				 raise ValueError('processmesh error: cannot work with 3D mesh in lat-lon coords')
			#we modify the mesh temporarily to a 2D mesh from which the 3D mesh was extruded
			 x=x2d
			 y=y2d
			 z=zeros(size(x2d))
			 elements=elements2d
	
	else:
		#Process mesh for plotting 
		if md.mesh.dimension()==2:
			is2d=1
		else:
			# process polycollection here for 3D plot
			is2d=0
	
	#units
	if options.exist('unit'):
		unit=options.getfieldvalue('unit')
		x=x*unit
		y=y*unit
		z=z*unit

	#is model a member of planet class? (workaround until planet class defined)
	if md.__class__.__name__!='model':
		isplanet=1
	else:
		isplanet=0

	return x,y,z,elements,is2d,isplanet
