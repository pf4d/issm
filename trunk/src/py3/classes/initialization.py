import numpy
from project3d import project3d
from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData
import MatlabFuncs as m

class initialization(object):
	"""
	INITIALIZATION class definition
	
	Usage:
	initialization=initialization();
	"""

	def __init__(self): # {{{
					
		self.vx            = float('NaN')
		self.vy            = float('NaN')
		self.vz            = float('NaN')
		self.vel           = float('NaN')
		self.pressure      = float('NaN')
		self.temperature   = float('NaN')
		self.waterfraction = float('NaN')
		self.watercolumn   = float('NaN')
		self.sediment_head = float('NaN')
		self.epl_head      = float('NaN')
		self.epl_thickness = float('NaN')

		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self): # {{{
		string='   initial field values:'
		string="%s\n%s"%(string,fielddisplay(self,'vx','x component of velocity [m/yr]'))
		string="%s\n%s"%(string,fielddisplay(self,'vy','y component of velocity [m/yr]'))
		string="%s\n%s"%(string,fielddisplay(self,'vz','z component of velocity [m/yr]'))
		string="%s\n%s"%(string,fielddisplay(self,'vel','velocity norm [m/yr]'))
		string="%s\n%s"%(string,fielddisplay(self,'pressure','pressure [Pa]'))
		string="%s\n%s"%(string,fielddisplay(self,'temperature','temperature [K]'))
		string="%s\n%s"%(string,fielddisplay(self,'waterfraction','fraction of water in the ice'))
		string="%s\n%s"%(string,fielddisplay(self,'watercolumn','thickness of subglacial water [m]'))
		string="%s\n%s"%(string,fielddisplay(self,'sediment_head','sediment water head of subglacial system [m]'))
		string="%s\n%s"%(string,fielddisplay(self,'epl_head','epl water head of subglacial system [m]'))
		string="%s\n%s"%(string,fielddisplay(self,'epl_thickness','thickness of the epl [m]'))

		return string
		#}}}
	def extrude(self,md): # {{{
		self.vx=project3d(md,'vector',self.vx,'type','node')
		self.vy=project3d(md,'vector',self.vy,'type','node')
		self.vz=project3d(md,'vector',self.vz,'type','node')
		self.vel=project3d(md,'vector',self.vel,'type','node')
		self.temperature=project3d(md,'vector',self.temperature,'type','node')
		self.waterfraction=project3d(md,'vector',self.waterfraction,'type','node')
		self.watercolumn=project3d(md,'vector',self.watercolumn,'type','node')
		self.sediment_head=project3d(md,'vector',self.sediment_head,'type','node','layer',1)
		self.epl_head=project3d(md,'vector',self.epl_head,'type','node','layer',1)
		self.epl_thickness=project3d(md,'vector',self.epl_thickness,'type','node','layer',1)

		#Lithostatic pressure by default
		self.pressure=md.constants.g*md.materials.rho_ice*(md.geometry.surface-md.mesh.z.reshape(-1,1))
		return self
	#}}}
	def setdefaultparameters(self): # {{{
		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{
		if StressbalanceAnalysisEnum() in analyses:
			if not numpy.any(numpy.logical_or(numpy.isnan(md.initialization.vx),numpy.isnan(md.initialization.vy))):
				md = checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
				md = checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
		if MasstransportAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
			md = checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
		if BalancethicknessAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
			md = checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
			#Triangle with zero velocity
			if numpy.any(numpy.logical_and(numpy.sum(numpy.abs(md.initialization.vx[md.mesh.elements-1]),axis=1)==0,\
			                               numpy.sum(numpy.abs(md.initialization.vy[md.mesh.elements-1]),axis=1)==0)):
				md.checkmessage("at least one triangle has all its vertices with a zero velocity")
		if ThermalAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
			md = checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
			md = checkfield(md,'fieldname','initialization.temperature','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
			if md.mesh.dimension()==3:
				md = checkfield(md,'fieldname','initialization.vz','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
			md = checkfield(md,'fieldname','initialization.pressure','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
			if (EnthalpyAnalysisEnum() in analyses and md.thermal.isenthalpy):
				md = checkfield(md,'fieldname','initialization.waterfraction','>=',0,'size',[md.mesh.numberofvertices])
				md = checkfield(md,'fieldname','initialization.watercolumn'  ,'>=',0,'size',[md.mesh.numberofvertices])
		if HydrologyShreveAnalysisEnum() in analyses:
			if hasattr(md.hydrology,'hydrologyshreve'):
				md = checkfield(md,'fieldname','initialization.watercolumn','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
		if HydrologyDCInefficientAnalysisEnum() in analyses:
			if hasattr(md.hydrology,'hydrologydc'):
				md = checkfield(md,'fieldname','initialization.sediment_head','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices,1])
		if HydrologyDCEfficientAnalysisEnum() in analyses:
			if hasattr(md.hydrology,'hydrologydc'):
				if md.hydrology.isefficientlayer==1:
					md = checkfield(md,'fieldname','initialization.epl_head','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices,1])
					md = checkfield(md,'fieldname','initialization.epl_thickness','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices,1])

		return md
	# }}}
	def marshall(self,md,fid):    # {{{

		yts=365.0*24.0*3600.0

		WriteData(fid,'data',self.vx,'format','DoubleMat','mattype',1,'enum',VxEnum(),'scale',1./yts)
		WriteData(fid,'data',self.vy,'format','DoubleMat','mattype',1,'enum',VyEnum(),'scale',1./yts)
		WriteData(fid,'data',self.vz,'format','DoubleMat','mattype',1,'enum',VzEnum(),'scale',1./yts)
		WriteData(fid,'data',self.pressure,'format','DoubleMat','mattype',1,'enum',PressureEnum())
		WriteData(fid,'data',self.temperature,'format','DoubleMat','mattype',1,'enum',TemperatureEnum())
		WriteData(fid,'data',self.waterfraction,'format','DoubleMat','mattype',1,'enum',WaterfractionEnum())
		WriteData(fid,'data',self.watercolumn,'format','DoubleMat','mattype',1,'enum',WatercolumnEnum())
		WriteData(fid,'data',self.sediment_head,'format','DoubleMat','mattype',1,'enum',SedimentHeadEnum())
		WriteData(fid,'data',self.epl_head,'format','DoubleMat','mattype',1,'enum',EplHeadEnum())
		WriteData(fid,'data',self.epl_thickness,'format','DoubleMat','mattype',1,'enum',HydrologydcEplThicknessEnum())

		
		if md.thermal.isenthalpy:
			tpmp = md.materials.meltingpoint - md.materials.beta*md.initialization.pressure;
			pos  = numpy.nonzero(md.initialization.temperature > tpmp)[0]
			enthalpy      = md.materials.heatcapacity*(md.initialization.temperature-md.constants.referencetemperature);
			enthalpy[pos] = md.materials.heatcapacity*tpmp[pos].reshape(-1,1) - md.constants.referencetemperature + md.materials.latentheat*md.initialization.waterfraction[pos].reshape(-1,1)
			WriteData(fid,'data',enthalpy,'format','DoubleMat','mattype',1,'enum',EnthalpyEnum());

	# }}}
