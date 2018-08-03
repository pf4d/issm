from fielddisplay import fielddisplay
from project3d import project3d
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData
import numpy

class basalforcings(object):
	"""
	BASAL FORCINGS class definition

	   Usage:
	      basalforcings=basalforcings();
	"""

	def __init__(self): # {{{
		self.groundedice_melting_rate  = float('NaN')
		self.floatingice_melting_rate  = float('NaN')
		self.geothermalflux            = float('NaN')

		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self): # {{{
		string="   basal forcings parameters:"

		string="%s\n%s"%(string,fielddisplay(self,"groundedice_melting_rate","basal melting rate (positive if melting) [m/yr]"))
		string="%s\n%s"%(string,fielddisplay(self,"floatingice_melting_rate","basal melting rate (positive if melting) [m/yr]"))
		string="%s\n%s"%(string,fielddisplay(self,"geothermalflux","geothermal heat flux [W/m^2]"))
		return string
		#}}}
	def extrude(self,md): # {{{
		self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1)
		self.floatingice_melting_rate=project3d(md,'vector',self.floatingice_melting_rate,'type','node','layer',1)
		self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','node','layer',1)    #bedrock only gets geothermal flux
		return self
	#}}}
	def initialize(self,md): # {{{

		if numpy.all(numpy.isnan(self.groundedice_melting_rate)):
			self.groundedice_melting_rate=numpy.zeros((md.mesh.numberofvertices,1))
			print("      no basalforcings.groundedice_melting_rate specified: values set as zero")

		if numpy.all(numpy.isnan(self.floatingice_melting_rate)):
			self.floatingice_melting_rate=numpy.zeros((md.mesh.numberofvertices,1))
			print("      no basalforcings.floatingice_melting_rate specified: values set as zero")

		return self
	#}}}
	def setdefaultparameters(self): # {{{
		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		if MasstransportAnalysisEnum() in analyses and not (solution==TransientSolutionEnum() and not md.transient.ismasstransport):
			md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1)
			md = checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'Inf',1,'timeseries',1)

		if BalancethicknessAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
			md = checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])

		if ThermalAnalysisEnum() in analyses and not (solution==TransientSolutionEnum() and not md.transient.isthermal):
			md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1)
			md = checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'Inf',1,'timeseries',1)
			md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'timeseries',1,'>=',0)

		return md
	# }}}
	def marshall(self,md,fid):    # {{{

		yts=365.0*24.0*3600.0

		WriteData(fid,'enum',BasalforcingsEnum(),'data',FloatingMeltRateEnum(),'format','Integer');
		WriteData(fid,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'fieldname','floatingice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'fieldname','geothermalflux','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1)
	# }}}
