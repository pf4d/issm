import numpy
from project3d import project3d
from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData
import MatlabFuncs as m

class thermal(object):
	"""
	THERMAL class definition

	   Usage:
	      thermal=thermal();
	"""

	def __init__(self): # {{{
		self.spctemperature    = float('NaN')
		self.penalty_threshold = 0
		self.stabilization     = 0
		self.reltol            = 0
		self.maxiter           = 0
		self.penalty_lock      = 0
		self.penalty_factor    = 0
		self.isenthalpy        = 0
		self.isdynamicbasalspc = 0;
		self.requested_outputs = []

		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self): # {{{
		string='   Thermal solution parameters:'
		string="%s\n%s"%(string,fielddisplay(self,'spctemperature','temperature constraints (NaN means no constraint) [K]'))
		string="%s\n%s"%(string,fielddisplay(self,'stabilization','0: no, 1: artificial_diffusivity, 2: SUPG'))
		string="%s\n%s"%(string,fielddisplay(self,'maxiter','maximum number of non linear iterations'))
		string="%s\n%s"%(string,fielddisplay(self,'reltol','relative tolerance criterion'))
		string="%s\n%s"%(string,fielddisplay(self,'penalty_lock','stabilize unstable thermal constraints that keep zigzagging after n iteration (default is 0, no stabilization)'))
		string="%s\n%s"%(string,fielddisplay(self,'penalty_threshold','threshold to declare convergence of thermal solution (default is 0)'))
		string="%s\n%s"%(string,fielddisplay(self,'isenthalpy','use an enthalpy formulation to include temperate ice (default is 0)'))
		string="%s\n%s"%(string,fielddisplay(self,'isdynamicbasalspc','enable dynamic setting of basal forcing. required for enthalpy formulation (default is 0)'))
		string="%s\n%s"%(string,fielddisplay(self,'requested_outputs','additional outputs requested'))
		return string
		#}}}
	def extrude(self,md): # {{{
		self.spctemperature=project3d(md,'vector',self.spctemperature,'type','node','layer',md.mesh.numberoflayers,'padding',numpy.nan)
		if isinstance(md.initialization.temperature,numpy.ndarray) and numpy.size(md.initialization.temperature,axis=0)==md.mesh.numberofvertices:
			self.spctemperature=numpy.nan*numpy.ones((md.mesh.numberofvertices,1))
			pos=numpy.nonzero(md.mesh.vertexonsurface)[0]
			self.spctemperature[pos]=md.initialization.temperature[pos]    #impose observed temperature on surface
		return self
	#}}}
	def defaultoutputs(self,md): # {{{

		if self.isenthalpy:
			return ['Enthalpy','Temperature','Waterfraction','Watercolumn','BasalforcingsGroundediceMeltingRate']
		else:
			return ['Temperature','BasalforcingsGroundediceMeltingRate']

	#}}}
	def setdefaultparameters(self): # {{{
		
		#Number of unstable constraints acceptable
		self.penalty_threshold=0

		#Type of stabilization used
		self.stabilization=1

		#Relative tolerance for the enthalpy convergence
		self.reltol=0.01

		#Maximum number of iterations
		self.maxiter=100

		#factor used to compute the values of the penalties: kappa=max(stiffness matrix)*10^penalty_factor
		self.penalty_factor=3

		#Should we use cold ice (default) or enthalpy formulation
		self.isenthalpy=0

		#will basal boundary conditions be set dynamically
		self.isdynamicbasalspc=0;

		#default output
		self.requested_outputs=['default']
		return self

	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		#Early return
		if (ThermalAnalysisEnum() not in analyses and EnthalpyAnalysisEnum() not in analyses) or (solution==TransientSolutionEnum() and not md.transient.isthermal):
			return md

		md = checkfield(md,'fieldname','thermal.stabilization','numel',[1],'values',[0,1,2])
		md = checkfield(md,'fieldname','thermal.spctemperature','Inf',1,'timeseries',1)
		if EnthalpyAnalysisEnum() in analyses and md.thermal.isenthalpy and md.mesh.dimension()==3:
			pos=numpy.nonzero(numpy.logical_not(numpy.isnan(md.thermal.spctemperature[0:md.mesh.numberofvertices])))
			replicate=numpy.tile(md.geometry.surface-md.mesh.z,(1,numpy.size(md.thermal.spctemperature,axis=1)))
			md = checkfield(md,'fieldname','thermal.spctemperature[numpy.nonzero(numpy.logical_not(numpy.isnan(md.thermal.spctemperature[0:md.mesh.numberofvertices,:])))]','<',md.materials.meltingpoint-md.materials.beta*md.materials.rho_ice*md.constants.g*replicate[pos],'message',"spctemperature should be below the adjusted melting point")
			md = checkfield(md,'fieldname','thermal.isenthalpy','numel',[1],'values',[0,1])
			md = checkfield(md,'fieldname','thermal.isdynamicbasalspc','numel',[1],'values',[0,1]);
			if(md.thermal.isenthalpy):
				if numpy.isnan(md.stressbalance.reltol):
					md.checkmessage("for a steadystate computation, thermal.reltol (relative convergence criterion) must be defined!")
				md = checkfield(md,'fieldname','thermal.reltol','>',0.,'message',"reltol must be larger than zero");
		md = checkfield(md,'fieldname','thermal.requested_outputs','stringrow',1)

		return md
	# }}}
	def marshall(self,md,fid):    # {{{
		WriteData(fid,'object',self,'fieldname','spctemperature','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'fieldname','penalty_threshold','format','Integer')
		WriteData(fid,'object',self,'fieldname','stabilization','format','Integer')
		WriteData(fid,'object',self,'fieldname','reltol','format','Double');
		WriteData(fid,'object',self,'fieldname','maxiter','format','Integer')
		WriteData(fid,'object',self,'fieldname','penalty_lock','format','Integer')
		WriteData(fid,'object',self,'fieldname','penalty_factor','format','Double')
		WriteData(fid,'object',self,'fieldname','isenthalpy','format','Boolean')
		WriteData(fid,'object',self,'fieldname','isdynamicbasalspc','format','Boolean');

		#process requested outputs
		outputs = self.requested_outputs
		indices = [i for i, x in enumerate(outputs) if x == 'default']
		if len(indices) > 0:
			outputscopy=outputs[0:max(0,indices[0]-1)]+self.defaultoutputs(md)+outputs[indices[0]+1:]
			outputs    =outputscopy
		WriteData(fid,'data',outputs,'enum',ThermalRequestedOutputsEnum(),'format','StringArray')
	# }}}
