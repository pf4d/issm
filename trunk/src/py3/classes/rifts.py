import numpy
from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData
from isnans import isnans
import MatlabFuncs as m

class rifts(object):
	"""
	RIFTS class definition

	   Usage:
	      rifts=rifts();
	"""

	def __init__(self): # {{{
		self.riftstruct     = []
		self.riftproperties = []

		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self): # {{{
		string='   rifts parameters:'

		string="%s\n%s"%(string,fielddisplay(self,'riftstruct','structure containing all rift information (vertices coordinates, segments, type of melange, ...)'))
		string="%s\n%s"%(string,fielddisplay(self,'riftproperties',''))
		return string
		#}}}
	def setdefaultparameters(self): # {{{
		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{
		if (not self.riftstruct) or numpy.any(isnans(self.riftstruct)):
			numrifts=0
		else:
			numrifts=len(self.riftstruct)

		if numrifts:
			if not m.strcmp(md.mesh.domaintype(),'2Dhorizontal'):
				md.checkmessage("models with rifts are only supported in 2d for now!")
			if not isinstance(self.riftstruct,list):
				md.checkmessage("rifts.riftstruct should be a structure!")
			if numpy.any(md.mesh.segmentmarkers>=2):
				#We have segments with rift markers, but no rift structure!
				md.checkmessage("model should be processed for rifts (run meshprocessrifts)!")
			for i,rift in enumerate(self.riftstruct):
				md = checkfield(md,'fieldname',"rifts.riftstruct[%d]['fill']" % i,'values',[WaterEnum(),AirEnum(),IceEnum(),MelangeEnum()])
		else:
			if self.riftstruct and numpy.any(numpy.logical_not(isnans(self.riftstruct))):
				md.checkmessage("riftstruct should be NaN since numrifts is 0!")

		return md
	# }}}
	def marshall(self,md,fid):    # {{{

		#Process rift info
		if (not self.riftstruct) or numpy.any(isnans(self.riftstruct)):
			numrifts=0
		else:
			numrifts=len(self.riftstruct)

		numpairs=0
		for rift in self.riftstruct:
			numpairs+=numpy.size(rift['penaltypairs'],axis=0)

		# 2 for nodes + 2 for elements+ 2 for  normals + 1 for length + 1 for fill + 1 for friction + 1 for fraction + 1 for fractionincrement + 1 for state.
		data=numpy.zeros((numpairs,12))
		count=0
		for rift in self.riftstruct:
			numpairsforthisrift=numpy.size(rift['penaltypairs'],0)
			data[count:count+numpairsforthisrift,0:7]=rift['penaltypairs']
			data[count:count+numpairsforthisrift,7]=rift['fill']
			data[count:count+numpairsforthisrift,8]=rift['friction']
			data[count:count+numpairsforthisrift,9]=rift['fraction']
			data[count:count+numpairsforthisrift,10]=rift['fractionincrement']
			data[count:count+numpairsforthisrift,11]=rift['state'].reshape(-1)
			count+=numpairsforthisrift

		WriteData(fid,'data',numrifts,'enum',RiftsNumriftsEnum(),'format','Integer')
		WriteData(fid,'data',data,'enum',RiftsRiftstructEnum(),'format','DoubleMat','mattype',3)
	# }}}
