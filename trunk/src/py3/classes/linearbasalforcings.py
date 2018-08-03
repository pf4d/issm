from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData
import numpy

class linearbasalforcings(object):
	"""
	LINEAR BASAL FORCINGS class definition

	   Usage:
	      basalforcings=linearbasalforcings();
	"""

	def __init__(self,*args): # {{{

		if not len(args):
			print('empty init')
			self.groundedice_melting_rate  = float('NaN')
			self.deepwater_melting_rate    = 0.
			self.deepwater_elevation       = 0.
			self.upperwater_elevation      = 0.
			self.geothermalflux            = float('NaN')

			#set defaults
			self.setdefaultparameters()
		elif len(args)==1 and args[0].__module__=='basalforcings':
			print('converting basalforings to linearbasalforcings')
			inv=args[0]
			self.groundedice_melting_rate  = inv.groundedice_melting_rate
			self.geothermalflux            = inv.geothermalflux
			self.deepwater_melting_rate    = 0.
			self.deepwater_elevation       = 0.
			self.upperwater_elevation      = 0.

			#set defaults
			self.setdefaultparameters()
		else:
			raise Exception('constructor not supported')

		#}}}
	def __repr__(self): # {{{
		string="   linear basal forcings parameters:"

		string="%s\n%s"%(string,fielddisplay(self,"groundedice_melting_rate","basal melting rate (positive if melting) [m/yr]"))
		string="%s\n%s"%(string,fielddisplay(self,"deepwater_melting_rate","basal melting rate (positive if melting applied for floating ice whith base < deepwater_elevation) [m/yr]"))
		string="%s\n%s"%(string,fielddisplay(self,"deepwater_elevation","elevation of ocean deepwater [m]"))
		string="%s\n%s"%(string,fielddisplay(self,"upperwater_elevation","elevation of ocean upper water [m]"))
		string="%s\n%s"%(string,fielddisplay(self,"geothermalflux","geothermal heat flux [W/m^2]"))
		return string
		#}}}
	def initialize(self,md): # {{{

		if numpy.all(numpy.isnan(self.groundedice_melting_rate)):
			self.groundedice_melting_rate=numpy.zeros((md.mesh.numberofvertices,1))
			print("      no basalforcings.groundedice_melting_rate specified: values set as zero")

		return self
	#}}}
	def setdefaultparameters(self): # {{{

		self.deepwater_melting_rate   = 50.0
		self.deepwater_elevation      = -800.0
		self.upperwater_elevation     = -400.0

		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		if MasstransportAnalysisEnum() in analyses and not (solution==TransientSolutionEnum() and not md.transient.ismasstransport):
			md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1)
			md = checkfield(md,'fieldname','basalforcings.deepwater_melting_rate','>=',0);
			md = checkfield(md,'fieldname','basalforcings.deepwater_elevation','<',md.basalforcings.upperwater_elevation);
			md = checkfield(md,'fieldname','basalforcings.upperwater_elevation','<',0);

		if BalancethicknessAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
			md = checkfield(md,'fieldname','basalforcings.deepwater_melting_rate','>=',0);
			md = checkfield(md,'fieldname','basalforcings.deepwater_elevation','<',md.basalforcings.upperwater_elevation);
			md = checkfield(md,'fieldname','basalforcings.upperwater_elevation','<',0);

		if ThermalAnalysisEnum() in analyses and not (solution==TransientSolutionEnum() and not md.transient.isthermal):
			md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1)
			md = checkfield(md,'fieldname','basalforcings.deepwater_melting_rate','>=',0);
			md = checkfield(md,'fieldname','basalforcings.deepwater_elevation','<',md.basalforcings.upperwater_elevation);
			md = checkfield(md,'fieldname','basalforcings.upperwater_elevation','<',0);
			md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'timeseries',1,'>=',0)

		return md
	# }}}
	def marshall(self,md,fid):    # {{{

		yts=365.0*24.0*3600.0

		floatingice_melting_rate = numpy.zeros((md.mesh.numberofvertices,1))
		pos=numpy.nonzero(md.geometry.base<=md.basalforcings.deepwater_elevation)
		floatingice_melting_rate[pos]=md.basalforcings.deepwater_melting_rate
		pos=numpy.nonzero(numpy.logical_and(md.geometry.base>md.basalforcings.deepwater_elevation,md.geometry.base<md.basalforcings.upperwater_elevation))
		floatingice_melting_rate[pos]=md.basalforcings.deepwater_melting_rate*(md.geometry.base[pos]-md.basalforcings.upperwater_elevation)/(md.basalforcings.deepwater_elevation-md.basalforcings.upperwater_elevation)

		WriteData(fid,'enum',BasalforcingsEnum(),'data',LinearFloatingMeltRateEnum(),'format','Integer');
		WriteData(fid,'object',self,'fieldname','groundedice_melting_rate','enum',BasalforcingsGroundediceMeltingRateEnum(),'format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'data',floatingice_melting_rate,'enum',BasalforcingsFloatingiceMeltingRateEnum(),'format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'fieldname','geothermalflux','enum',BasalforcingsGeothermalfluxEnum(),'format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'fieldname','deepwater_melting_rate','enum',BasalforcingsDeepwaterMeltingRateEnum(),'format','Double','scale',1./yts)
		WriteData(fid,'object',self,'fieldname','deepwater_elevation','enum',BasalforcingsDeepwaterElevationEnum(),'format','Double')
		WriteData(fid,'object',self,'fieldname','upperwater_elevation','enum',BasalforcingsUpperwaterElevationEnum(),'format','Double')
	# }}}
