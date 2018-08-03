import os.path
import numpy
import copy
from NodeConnectivity import NodeConnectivity
from ElementConnectivity import ElementConnectivity
from mesh2d import mesh2d
from mesh3dprisms import mesh3dprisms
import MatlabFuncs as m

def contourenvelope(md,*args):
	"""
	CONTOURENVELOPE - build a set of segments enveloping a contour .exp

	   Usage:
	      segments=contourenvelope(md,varargin)

	   Example:
	      segments=contourenvelope(md,'Stream.exp');
	      segments=contourenvelope(md);
	"""

	#some checks
	if len(args)>1:
		raise RuntimeError("contourenvelope error message: bad usage")

	if len(args)==1:
		flags=args[0]

		if   isinstance(flags,str):
			file=flags
			if not os.path.exists(file):
				raise IOError("contourenvelope error message: file '%s' not found" % file)
			isfile=1
		elif isinstance(flags,(bool,int,float)):
			#do nothing for now
			isfile=0
		else:
			raise TypeError("contourenvelope error message:  second argument should be a file or an elements flag")

	#Now, build the connectivity tables for this mesh.
	#Computing connectivity
	if numpy.size(md.mesh.vertexconnectivity,axis=0)!=md.mesh.numberofvertices and numpy.size(md.mesh.vertexconnectivity,axis=0)!=md.mesh.numberofvertices2d:
		[md.mesh.vertexconnectivity]=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices)
	if numpy.size(md.mesh.elementconnectivity,axis=0)!=md.mesh.numberofelements and numpy.size(md.mesh.elementconnectivity,axis=0)!=md.mesh.numberofelements2d:
		[md.mesh.elementconnectivity]=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity)

	#get nodes inside profile
	elementconnectivity=copy.deepcopy(md.mesh.elementconnectivity)
	if md.mesh.dimension()==2:
		elements=copy.deepcopy(md.mesh.elements)
		x=copy.deepcopy(md.mesh.x)
		y=copy.deepcopy(md.mesh.y)
		numberofvertices=copy.deepcopy(md.mesh.numberofvertices)
		numberofelements=copy.deepcopy(md.mesh.numberofelements)
	else:
		elements=copy.deepcopy(md.mesh.elements2d)
		x=copy.deepcopy(md.mesh.x2d)
		y=copy.deepcopy(md.mesh.y2d)
		numberofvertices=copy.deepcopy(md.mesh.numberofvertices2d)
		numberofelements=copy.deepcopy(md.mesh.numberofelements2d)

	if len(args)==1:

		if isfile:
			#get flag list of elements and nodes inside the contour
			nodein=ContourToMesh(elements,x,y,file,'node',1)
			elemin=(numpy.sum(nodein(elements),axis=1)==numpy.size(elements,axis=1))
			#modify element connectivity
			elemout=numpy.nonzero(numpy.logical_not(elemin))[0]
			elementconnectivity[elemout,:]=0
			elementconnectivity[numpy.nonzero(m.ismember(elementconnectivity,elemout+1))]=0
		else:
			#get flag list of elements and nodes inside the contour
			nodein=numpy.zeros(numberofvertices)
			elemin=numpy.zeros(numberofelements)

			pos=numpy.nonzero(flags)
			elemin[pos]=1
			nodein[elements[pos,:]-1]=1

			#modify element connectivity
			elemout=numpy.nonzero(numpy.logical_not(elemin))[0]
			elementconnectivity[elemout,:]=0
			elementconnectivity[numpy.nonzero(m.ismember(elementconnectivity,elemout+1))]=0

	#Find element on boundary
	#First: find elements on the boundary of the domain
	flag=copy.deepcopy(elementconnectivity)
	if len(args)==1:
		flag[numpy.nonzero(flag)]=elemin[flag[numpy.nonzero(flag)]]
	elementonboundary=numpy.logical_and(numpy.prod(flag,axis=1)==0,numpy.sum(flag,axis=1)>0)

	#Find segments on boundary
	pos=numpy.nonzero(elementonboundary)[0]
	num_segments=numpy.size(pos)
	segments=numpy.zeros((num_segments*3,3),int)
	count=0

	for el1 in pos:
		els2=elementconnectivity[el1,numpy.nonzero(elementconnectivity[el1,:])[0]]-1
		if numpy.size(els2)>1:
			flag=numpy.intersect1d(numpy.intersect1d(elements[els2[0],:],elements[els2[1],:]),elements[el1,:])
			nods1=elements[el1,:]
			nods1=numpy.delete(nods1,numpy.nonzero(nods1==flag))
			segments[count,:]=[nods1[0],nods1[1],el1+1]

			ord1=numpy.nonzero(nods1[0]==elements[el1,:])[0][0]
			ord2=numpy.nonzero(nods1[1]==elements[el1,:])[0][0]

			#swap segment nodes if necessary
			if ( (ord1==0 and ord2==1) or (ord1==1 and ord2==2) or (ord1==2 and ord2==0) ):
				temp=segments[count,0]
				segments[count,0]=segments[count,1]
				segments[count,1]=temp
			segments[count,0:2]=numpy.flipud(segments[count,0:2])
			count+=1
		else:
			nods1=elements[el1,:]
			flag=numpy.setdiff1d(nods1,elements[els2,:])
			for j in range(0,3):
				nods=numpy.delete(nods1,j)
				if numpy.any(m.ismember(flag,nods)):
					segments[count,:]=[nods[0],nods[1],el1+1]
					ord1=numpy.nonzero(nods[0]==elements[el1,:])[0][0]
					ord2=numpy.nonzero(nods[1]==elements[el1,:])[0][0]
					if ( (ord1==0 and ord2==1) or (ord1==1 and ord2==2) or (ord1==2 and ord2==0) ):
						temp=segments[count,0]
						segments[count,0]=segments[count,1]
						segments[count,1]=temp
					segments[count,0:2]=numpy.flipud(segments[count,0:2])
					count+=1
	segments=segments[0:count,:]

	return segments

