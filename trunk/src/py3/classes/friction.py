from fielddisplay import fielddisplay
from project3d import project3d
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData

class friction(object):
	"""
	FRICTION class definition

	   Usage:
	      friction=friction()
	"""

	def __init__(self): # {{{
		self.coefficient = float('NaN')
		self.p           = float('NaN')
		self.q           = float('NaN')

		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self): # {{{
		string="Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s-1) * u_b,\n(effective stress Neff=rho_ice*g*thickness+rho_water*g*bed, r=q/p and s=1/p)"

		string="%s\n%s"%(string,fielddisplay(self,"coefficient","friction coefficient [SI]"))
		string="%s\n%s"%(string,fielddisplay(self,"p","p exponent"))
		string="%s\n%s"%(string,fielddisplay(self,"q","q exponent"))
		return string
		#}}}
	def extrude(self,md): # {{{
		self.coefficient=project3d(md,'vector',self.coefficient,'type','node','layer',1)
		self.p=project3d(md,'vector',self.p,'type','element')
		self.q=project3d(md,'vector',self.q,'type','element')
		return self
	#}}}
	def setdefaultparameters(self): # {{{
		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		#Early return
		if StressbalanceAnalysisEnum() not in analyses and ThermalAnalysisEnum() not in analyses:
			return md

		md = checkfield(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1,'Inf',1)
		md = checkfield(md,'fieldname','friction.q','NaN',1,'Inf',1,'size',[md.mesh.numberofelements])
		md = checkfield(md,'fieldname','friction.p','NaN',1,'Inf',1,'size',[md.mesh.numberofelements])

		return md
	# }}}
	def marshall(self,md,fid):    # {{{
		WriteData(fid,'enum',FrictionLawEnum(),'data',1,'format','Integer')
		WriteData(fid,'object',self,'fieldname','coefficient','format','DoubleMat','mattype',1)
		WriteData(fid,'object',self,'fieldname','p','format','DoubleMat','mattype',2)
		WriteData(fid,'object',self,'fieldname','q','format','DoubleMat','mattype',2)
	# }}}
