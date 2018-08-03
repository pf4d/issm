from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData

class hydrologyshreve(object):
	"""
	HYDROLOGYSHREVE class definition

	   Usage:
	      hydrologyshreve=hydrologyshreve();
	"""

	def __init__(self): # {{{
		self.spcwatercolumn = float('NaN')
		self.stabilization  = 0

		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self): # {{{
		
		string='   hydrologyshreve solution parameters:'
		string="%s\n%s"%(string,fielddisplay(self,'spcwatercolumn','water thickness constraints (NaN means no constraint) [m]'))
		string="%s\n%s"%(string,fielddisplay(self,'stabilization','artificial diffusivity (default is 1). can be more than 1 to increase diffusivity.'))
		return string
		#}}}
	def extrude(self,md): # {{{
		return self
	#}}}
	def setdefaultparameters(self): # {{{
		
		#Type of stabilization to use 0:nothing 1:artificial_diffusivity
		self.stabilization=1

		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{
		
		#Early return
		if HydrologyShreveAnalysisEnum() not in analyses:
			return md

		md = checkfield(md,'fieldname','hydrology.spcwatercolumn','Inf',1,'timeseries',1)
		md = checkfield(md,'fieldname','hydrology.stabilization','>=',0)

		return md
	# }}}
	def marshall(self,md,fid):    # {{{
		WriteData(fid,'enum',HydrologyModelEnum(),'data',HydrologyshreveEnum(),'format','Integer');
		WriteData(fid,'object',self,'fieldname','spcwatercolumn','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'fieldname','stabilization','format','Double')
	# }}}
