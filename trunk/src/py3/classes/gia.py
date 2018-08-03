from fielddisplay import fielddisplay
from project3d import project3d
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData

class gia(object):
	"""
	GIA class definition

	   Usage:
	      gia=gia();
	"""

	def __init__(self): # {{{
		self.mantle_viscosity              = float('NaN');
		self.lithosphere_thickness         = float('NaN');
		self.cross_section_shape           = 0;
	
		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self): # {{{
		
		string='   gia solution parameters:' 
		
		string="%s\n%s"%(string,fielddisplay(self,'mantle_viscosity','mantle viscosity constraints (NaN means no constraint) (Pa s)'))
		string="%s\n%s"%(string,fielddisplay(self,'lithosphere_thickness','lithosphere thickness constraints (NaN means no constraint) (m)'))
		string="%s\n%s"%(string,fielddisplay(self,'cross_section_shape',"1: square-edged, 2: elliptical-edged surface"))
		return string
		#}}}
	def extrude(self,md): # {{{
		self.mantle_viscosity=project3d(md,'vector',self.mantle_viscosity,'type','node')
		self.lithosphere_thickness=project3d(md,'vector',self.lithosphere_thickness,'type','node')
		return self
	#}}}
	def setdefaultparameters(self): # {{{

		self.cross_section_shape=1; 

		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		# Early return 
		if (GiaAnalysisEnum() not in  analyses):
			return md 
		
		md = checkfield(md,'fieldname','gia.mantle_viscosity','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices,1],'>',0)
		md = checkfield(md,'fieldname','gia.lithosphere_thickness','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices,1],'>',0)
		md = checkfield(md,'fieldname','gia.cross_section_shape','numel',[1],'values',[1,2])

		#be sure that if we are running a masstransport ice flow model coupled with gia, that thickness forcings 
		#are not provided into the future.

		return md
	# }}}
	def marshall(self,md,fid):    # {{{

		WriteData(fid,'object',self,'fieldname','mantle_viscosity','format','DoubleMat','mattype',1);
		WriteData(fid,'object',self,'fieldname','lithosphere_thickness','format','DoubleMat','mattype',1,'scale',10.**3.);
		WriteData(fid,'object',self,'fieldname','cross_section_shape','format','Integer');
	# }}}
