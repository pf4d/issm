from fielddisplay import fielddisplay
from project3d import project3d
from EnumDefinitions import *
from StringToEnum import StringToEnum
from checkfield import checkfield
from WriteData import WriteData

class calving(object):
	"""
	CALVING class definition

	   Usage:
	      calving=calving();
	"""

	def __init__(self): # {{{

		self.stabilization = 0
		self.spclevelset   = float('NaN')
		self.calvingrate   = float('NaN')
		self.meltingrate   = float('NaN')

		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self): # {{{
		string='   Calving parameters:'
		string="%s\n%s"%(string,fielddisplay(self,'stabilization','0: no, 1: artificial_diffusivity, 2: streamline upwinding'))
		string="%s\n%s"%(string,fielddisplay(self,'spclevelset','levelset constraints (NaN means no constraint)'))
		string="%s\n%s"%(string,fielddisplay(self,'calvingrate','calving rate at given location [m/a]'))
		string="%s\n%s"%(string,fielddisplay(self,'meltingrate','melting rate at given location [m/a]'))

		return string
		#}}}
	def extrude(self,md): # {{{
		self.spclevelset=project3d(md,'vector',self.spclevelset,'type','node')
		self.calvingrate=project3d(md,'vector',self.calvingrate,'type','node')
		self.meltingrate=project3d(md,'vector',self.meltingrate,'type','node')
		return self
	#}}}
	def setdefaultparameters(self): # {{{

		#stabilization = 2 by default
		self.stabilization = 2

		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		#Early return
		if (solution!=TransientSolutionEnum()) or (not md.transient.iscalving):
			return md

		md = checkfield(md,'fieldname','calving.spclevelset','Inf',1,'timeseries',1)
		md = checkfield(md,'fieldname','calving.stabilization','values',[0,1,2]);
		md = checkfield(md,'fieldname','calving.calvingrate','>=',0,'timeseries',1,'NaN',1,'Inf',1);
		md = checkfield(md,'fieldname','calving.meltingrate','>=',0,'timeseries',1,'NaN',1,'Inf',1);

		return md
	# }}}
	def marshall(self,md,fid):    # {{{

		yts=365.*24.*3600.

		WriteData(fid,'enum',CalvingLawEnum(),'data',DefaultCalvingEnum(),'format','Integer');
		WriteData(fid,'enum',LevelsetStabilizationEnum(),'data',self.stabilization,'format','Integer');
		WriteData(fid,'enum',SpclevelsetEnum(),'data',self.spclevelset,'format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1);
		WriteData(fid,'object',self,'fieldname','calvingrate','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'scale',1./yts)
		WriteData(fid,'object',self,'fieldname','meltingrate','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'scale',1./yts)
	# }}}
