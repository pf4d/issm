import numpy
from GetNodalFunctionsCoeff import GetNodalFunctionsCoeff
from GetAreas import GetAreas
import MatlabFuncs as m

def ComputeHessian(index,x,y,field,type):
	"""
	COMPUTEHESSIAN - compute hessian matrix from a field

	   Compute the hessian matrix of a given field
	   return the three components Hxx Hxy Hyy
	   for each element or each node

	   Usage:
	      hessian=ComputeHessian(index,x,y,field,type)

	   Example:
	      hessian=ComputeHessian(md.mesh.elements,md.mesh.x,md.mesh.y,md.inversion.vel_obs,'node')
	"""

	#some variables
	numberofnodes=numpy.size(x)
	numberofelements=numpy.size(index,axis=0)

	#some checks
	if numpy.size(field)!=numberofnodes and numpy.size(field)!=numberofelements:
		raise TypeError("ComputeHessian error message: the given field size not supported yet")
	if not m.strcmpi(type,'node') and not m.strcmpi(type,'element'):
		raise TypeError("ComputeHessian error message: only 'node' or 'element' type supported yet")

	#initialization
	line=index.reshape(-1,order='F')
	linesize=3*numberofelements

	#get areas and nodal functions coefficients N(x,y)=alpha x + beta y + gamma 
	[alpha,beta,dum]=GetNodalFunctionsCoeff(index,x,y)
	areas=GetAreas(index,x,y)

	#compute weights that hold the volume of all the element holding the node i
	weights=m.sparse(line,numpy.ones((linesize,1)),numpy.tile(areas.reshape(-1,1),(3,1)),numberofnodes,1)

	#compute field on nodes if on elements
	if numpy.size(field,axis=0)==numberofelements:
		field=m.sparse(line,numpy.ones((linesize,1)),numpy.tile(areas*field,(3,1)),numberofnodes,1)/weights

	#Compute gradient for each element
	grad_elx=numpy.sum(field[index-1,0]*alpha,axis=1) 
	grad_ely=numpy.sum(field[index-1,0]*beta,axis=1)

	#Compute gradient for each node (average of the elements around)
	gradx=m.sparse(line,numpy.ones((linesize,1)),numpy.tile((areas*grad_elx).reshape(-1,1),(3,1)),numberofnodes,1)
	grady=m.sparse(line,numpy.ones((linesize,1)),numpy.tile((areas*grad_ely).reshape(-1,1),(3,1)),numberofnodes,1)
	gradx=gradx/weights
	grady=grady/weights

	#Compute hessian for each element
	hessian=numpy.hstack((numpy.sum(gradx[index-1,0]*alpha,axis=1).reshape(-1,1),numpy.sum(grady[index-1,0]*alpha,axis=1).reshape(-1,1),numpy.sum(grady[index-1,0]*beta,axis=1).reshape(-1,1)))

	if m.strcmpi(type,'node'):
		#Compute Hessian on the nodes (average of the elements around)
		hessian=numpy.hstack((m.sparse(line,numpy.ones((linesize,1)),numpy.tile((areas*hessian[:,0]).reshape(-1,1),(3,1)),numberofnodes,1)/weights, \
			m.sparse(line,numpy.ones((linesize,1)),numpy.tile((areas*hessian[:,1]).reshape(-1,1),(3,1)),numberofnodes,1)/weights, \
			m.sparse(line,numpy.ones((linesize,1)),numpy.tile((areas*hessian[:,2]).reshape(-1,1),(3,1)),numberofnodes,1)/weights ))

	return hessian

