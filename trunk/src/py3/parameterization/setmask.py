import numpy
import os
from model import model
from FlagElements import FlagElements
from pairoptions import pairoptions
from ContourToMesh import ContourToMesh

def setmask(md, floatingicename, groundedicename, **kwargs):
	"""
	SETMASK - establish boundaries between grounded and floating ice.

	   By default, ice is considered grounded. The contour floatingicename defines nodes 
	   for which ice is floating. The contour groundedicename defines nodes inside an floatingice, 
	   that are grounded (ie: ice rises, islands, etc ...)
	   All input files are in the Argus format (extension .exp).

	   Usage:
	      md=setmask(md,floatingicename,groundedicename)

	   Examples:
	      md=setmask(md,'all','');
	      md=setmask(md,'Iceshelves.exp','Islands.exp');
	"""
	#some checks on list of arguments
	if not isinstance(md,model):
		raise TypeError("setmask error message")

	#process options
	options=pairoptions(**kwargs)

	#Get assigned fields
	x = md.mesh.x
	y = md.mesh.y
	elements = md.mesh.elements

	#Assign elementonfloatingice, elementongroundedice, vertexongroundedice and vertexonfloatingice. Only change at your own peril! This is synchronized heavily with the GroundingLineMigration module. {{{
	elementonfloatingice = FlagElements(md, floatingicename)
	elementongroundedice = FlagElements(md, groundedicename) 

	#Because groundedice nodes and elements can be included into an floatingice, we need to update. Remember, all the previous 
	#arrays come from domain outlines that can intersect one another: 

	elementonfloatingice = numpy.logical_and(elementonfloatingice,numpy.logical_not(elementongroundedice))
	elementongroundedice = numpy.logical_not(elementonfloatingice)

	#the order here is important. we choose vertexongroundedice as default on the grounding line.
	vertexonfloatingice = numpy.zeros(md.mesh.numberofvertices,'bool')
	vertexongroundedice = numpy.zeros(md.mesh.numberofvertices,'bool')
	vertexongroundedice[md.mesh.elements[numpy.nonzero(elementongroundedice),:]-1]=True
	vertexonfloatingice[numpy.nonzero(numpy.logical_not(vertexongroundedice))]=True
	#}}}

	#level sets
	md.mask.groundedice_levelset = -1.*numpy.ones(md.mesh.numberofvertices)
	md.mask.groundedice_levelset[md.mesh.elements[numpy.nonzero(elementongroundedice),:]-1]=1.

	if(len(kwargs)):
		md.mask.ice_levelset = 1.*numpy.ones(md.mesh.numberofvertices)
		icedomainfile = options.getfieldvalue('icedomain','none')
		if not os.path.exists(icedomainfile):
			raise IOError("setmask error message: ice domain file '%s' not found." % icedomainfile)
		#use contourtomesh to set ice values inside ice domain
		[vertexinsideicedomain,elementinsideicedomain]=ContourToMesh(elements,x,y,icedomainfile,'node',1)
		md.mask.ice_levelset[numpy.nonzero(vertexinsideicedomain)[0]] = -1.
	else:
		md.mask.ice_levelset = -1.*numpy.ones(md.mesh.numberofvertices)

	return md
