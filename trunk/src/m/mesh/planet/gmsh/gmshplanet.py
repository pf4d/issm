from MatlabFuncs import *
import numpy as np
from numpy import *
from pairoptions import *
from mesh3dsurface import *
import subprocess

def gmshplanet(*varargin):
#GMSHPLANET - mesh generation for a sphere. Very specific code for gmsh. From demo/sphere.geo
#
#   Available options (for more details see ISSM website http://issm.jpl.nasa.gov/):
#
#   - radius:             radius of the planet in km
#   - resolution:         resolution in km
#   - refine:             provide mesh
#   - refinemetric:       mesh quantity to specify resolution
#
#   Returns 'mesh3dsurface' type mesh
#
#   Examples:
#      md.mesh=gmshplanet('radius',6000,'resolution',100);
#      md.mesh=gmshplanet('radius',6000,'resolution',100);

	#process options
	options=pairoptions(*varargin)
	#options=deleteduplicates(options,1)

	#recover parameters:
	radius=options.getfieldvalue('radius')*1000
	resolution=options.getfieldvalue('resolution')*1000

	#initialize mesh: 
	mesh=mesh3dsurface()

	#create .geo file:  {{{
	fid=open('sphere.geo','w')

	fid.write('Mesh.Algorithm = 1;\n')
	if options.exist('refine'):
		fid.write('Mesh.Algorithm = 7;\n')
		fid.write('Mesh.CharacteristicLengthFromPoints= 0;\n')
		fid.write('Mesh.SmoothRatio= 3;\n')
		fid.write('Mesh.RemeshAlgorithm= 1;\n')
	fid.write('resolution=%g;\n'%resolution)
	fid.write('radius=%g;\n'%radius)
	fid.write('Point(1) = {0.0,0.0,0.0,resolution};\n')
	fid.write('Point(2) = {radius,0.0,0.0,resolution};\n')
	fid.write('Point(3) = {0,radius,0.0,resolution};\n')
	fid.write('Circle(1) = {2,1,3};\n')
	fid.write('Point(4) = {-radius,0,0.0,resolution};\n')
	fid.write('Point(5) = {0,-radius,0.0,resolution};\n')
	fid.write('Circle(2) = {3,1,4};\n')
	fid.write('Circle(3) = {4,1,5};\n')
	fid.write('Circle(4) = {5,1,2};\n')
	fid.write('Point(6) = {0,0,-radius,resolution};\n')
	fid.write('Point(7) = {0,0,radius,resolution};\n')
	fid.write('Circle(5) = {3,1,6};\n')
	fid.write('Circle(6) = {6,1,5};\n')
	fid.write('Circle(7) = {5,1,7};\n')
	fid.write('Circle(8) = {7,1,3};\n')
	fid.write('Circle(9) = {2,1,7};\n')
	fid.write('Circle(10) = {7,1,4};\n')
	fid.write('Circle(11) = {4,1,6};\n')
	fid.write('Circle(12) = {6,1,2};\n')
	fid.write('Line Loop(13) = {2,8,-10};\n')
	fid.write('Ruled Surface(14) = {13};\n')
	fid.write('Line Loop(15) = {10,3,7};\n')
	fid.write('Ruled Surface(16) = {15};\n')
	fid.write('Line Loop(17) = {-8,-9,1};\n')
	fid.write('Ruled Surface(18) = {17};\n')
	fid.write('Line Loop(19) = {-11,-2,5};\n')
	fid.write('Ruled Surface(20) = {19};\n')
	fid.write('Line Loop(21) = {-5,-12,-1};\n')
	fid.write('Ruled Surface(22) = {21};\n')
	fid.write('Line Loop(23) = {-3,11,6};\n')
	fid.write('Ruled Surface(24) = {23};\n')
	fid.write('Line Loop(25) = {-7,4,9};\n')
	fid.write('Ruled Surface(26) = {25};\n')
	fid.write('Line Loop(27) = {-4,12,-6};\n')
	fid.write('Ruled Surface(28) = {27};\n')
	fid.write('Surface Loop(29) = {28,26,16,14,20,24,22,18};\n')
	fid.write('Volume(30) = {29};\n')
	fid.write('Physical Surface(1) = {28,26,16,14,20,24,22,18};\n')
	fid.write('Physical Volume(2) = 30;\n')
	fid.close()
	#}}}

	if options.exist('refine'):
		meshini=options.getfieldvalue('refine')
		metric=options.getfieldvalue('refinemetric')

		#create .pos file with existing mesh and refining metric:  {{{
		fid=open('sphere.pos','w')

		fid.write('View "background mesh" [;\n')
		for i in range(meshini.numberofelements):
			fid.write('ST(%g,%g,%g,%g,%g,%g,%g,%g,%g)[%g,%g,%g];\n',
								meshini.x(meshini.elements(i,0)), meshini.y(meshini.elements(i,0)), meshini.z(meshini.elements(i,0)),
								meshini.x(meshini.elements(i,1)), meshini.y(meshini.elements(i,1)), meshini.z(meshini.elements(i,1)),
								meshini.x(meshini.elements(i,2)), meshini.y(meshini.elements(i,2)), meshini.z(meshini.elements(i,2)),
								metric(meshini.elements(i,0)), metric(meshini.elements(i,1)), metric(meshini.elements(i,2)))
		fid.write('];\n')
		
		fid.close()
		# }}}

	#call gmsh
	if options.exist('refine'):
		subprocess.call('gmsh -tol 1e-8 -2 sphere.geo -bgm sphere.pos',shell=True)
	else:
		#call gmsh
		subprocess.call('gmsh -tol 1e-8 -2 sphere.geo',shell=True)

	#import mesh:  {{{
	fid=open('sphere.msh','r')

	#Get Mesh format
	A=fid.readline().strip()
	if not strcmp(A,'$MeshFormat'):
		raise RuntimeError(['Expecting $MeshFormat (', A, ')'])

	A=fid.readline().split()
	A=fid.readline().strip()
	if not strcmp(A,'$EndMeshFormat'):
		raise RuntimeError(['Expecting $EndMeshFormat (', A, ')'])

	#Nodes
	A=fid.readline().strip()
	if not strcmp(A,'$Nodes'):
		raise RuntimeError(['Expecting $Nodes (', A, ')'])

	mesh.numberofvertices=int(fid.readline().strip())
	mesh.x=np.empty(mesh.numberofvertices)
	mesh.y=np.empty(mesh.numberofvertices)
	mesh.z=np.empty(mesh.numberofvertices)
	for i in range(mesh.numberofvertices):
		A=fid.readline().split()
		mesh.x[i]=float(A[1])
		mesh.y[i]=float(A[2])
		mesh.z[i]=float(A[3])

	A=fid.readline().strip()
	if not strcmp(A,'$EndNodes'):
		raise RuntimeError(['Expecting $EndNodes (', A, ')'])

	#Elements
	A=fid.readline().strip()
	if not strcmp(A,'$Elements'):
		raise RuntimeError(['Expecting $Elements (', A, ')'])
	mesh.numberofelements=int(fid.readline().strip())
	mesh.elements=np.zeros([mesh.numberofelements,3])
	for i in range(mesh.numberofelements):
		A=fid.readline().split()
		mesh.elements[i]=[int(A[5]),int(A[6]),int(A[7])]
	mesh.elements=mesh.elements.astype(int)
	A=fid.readline().strip()
	if not strcmp(A,'$EndElements'):
		raise RuntimeError(['Expecting $EndElements (', A, ')'])
	fid.close() 
	#}}}

	#figure out other fields in mesh3dsurface: 
	mesh.r=np.sqrt(mesh.x**2+mesh.y**2+mesh.z**2)
	mesh.lat=np.arcsin(mesh.z/mesh.r)/np.pi*180
	mesh.long=np.arctan2(mesh.y,mesh.x)/np.pi*180

	#erase files: 
	subprocess.call('rm -rf sphere.geo sphere.msh sphere.pos',shell=True)

	#return mesh: 
	return mesh
