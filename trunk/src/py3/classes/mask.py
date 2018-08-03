import numpy
from fielddisplay import fielddisplay
from project3d import project3d
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData
import MatlabFuncs as m

class mask(object):
	"""
	MASK class definition

	   Usage:
	      mask=mask();
	"""

	def __init__(self): # {{{
		self.ice_levelset         = float('NaN')
		self.groundedice_levelset = float('NaN')

		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self): # {{{
		string="   masks:"

		string="%s\n%s"%(string,fielddisplay(self,"groundedice_levelset","is ice grounded ? grounded ice if > 0, grounding line position if = 0, floating ice if < 0"))
		string="%s\n%s"%(string,fielddisplay(self,"ice_levelset","presence of ice if < 0, icefront position if = 0, no ice if > 0"))
		return string
		#}}}
	def extrude(self,md): # {{{
		self.ice_levelset=project3d(md,'vector',self.ice_levelset,'type','node')
		self.groundedice_levelset=project3d(md,'vector',self.groundedice_levelset,'type','node')
		return self
	#}}}
	def setdefaultparameters(self): # {{{
		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		md = checkfield(md,'fieldname','mask.ice_levelset'        ,'size',[md.mesh.numberofvertices])
		isice=numpy.array(md.mask.ice_levelset<=0,int)
		if numpy.sum(isice)==0:
			raise TypeError("no ice present in the domain")

		icefront=numpy.sum(md.mask.ice_levelset[md.mesh.elements-1]==0,axis=1)
		if (max(icefront)==3 and m.strcmp(md.mesh.elementtype(),'Tria')) or (max(icefront==6) and m.strcmp(md.mesh.elementtype(),'Penta')):
			raise TypeError("At least one element has all nodes on ice front, change md.mask.ice_levelset to fix it")

		return md
	# }}}
	def marshall(self,md,fid):    # {{{
		WriteData(fid,'object',self,'fieldname','groundedice_levelset','format','DoubleMat','mattype',1)
		WriteData(fid,'object',self,'fieldname','ice_levelset','format','DoubleMat','mattype',1)

		# get mask of vertices of elements with ice
		isice=numpy.array(md.mask.ice_levelset<0.,int)
		vlist = numpy.zeros((md.mesh.numberofvertices,1), dtype=int)
		pos=numpy.nonzero(numpy.sum(isice[md.mesh.elements-1],axis=1))[0]
		vlist[md.mesh.elements[pos,:]-1]=1
		WriteData(fid,'data',vlist,'enum',IceMaskNodeActivationEnum(),'format','DoubleMat','mattype',1);
	# }}}
