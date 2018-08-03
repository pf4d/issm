from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import *
from project3d import *
from WriteData import *

class SMBcomponents(object):
	"""
	SMBcomponents Class definition

	   Usage:
	      SMBcomponents=SMBcomponents();
	"""

	def __init__(self): # {{{
		self.accumulation = float('NaN')
		self.runoff = float('NaN')
		self.evaporation = float('NaN')
		self.requested_outputs      = []
		#}}}
	def __repr__(self): # {{{
		string="   surface forcings parameters (SMB=accumulation-runoff-evaporation) :"
		string="%s\n%s"%(string,fielddisplay(self,'accumulation','accumulated snow [m/yr ice eq]'))
		string="%s\n%s"%(string,fielddisplay(self,'runoff','amount of ice melt lost from the ice column [m/yr ice eq]'))
		string="%s\n%s"%(string,fielddisplay(self,'evaporation','mount of ice lost to evaporative processes [m/yr ice eq]'))
		string="%s\n%s"%(string,fielddisplay(self,'requested_outputs','additional outputs requested'))
		return string
		#}}}
	def extrude(self,md): # {{{

		self.mass_balance=project3d(md,'vector',self.accumulation,'type','node');
		self.mass_balance=project3d(md,'vector',self.runoff,'type','node');
		self.mass_balance=project3d(md,'vector',self.evaporation,'type','node');
		return self
	#}}}
	def defaultoutputs(self,md): # {{{
		return []
	#}}}
	def initialize(self,md): # {{{

		if numpy.all(numpy.isnan(self.accumulation)):
			self.accumulation=numpy.zeros((md.mesh.numberofvertices,1))
			print("      no SMB.accumulation specified: values set as zero")

		if numpy.all(numpy.isnan(self.runoff)):
			self.runoff=numpy.zeros((md.mesh.numberofvertices,1))
			print("      no SMB.runoff specified: values set as zero")

		if numpy.all(numpy.isnan(self.evaporation)):
			self.evaporation=numpy.zeros((md.mesh.numberofvertices,1))
			print("      no SMB.evaporation specified: values set as zero")

		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		if MasstransportAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','smb.accumulation','timeseries',1,'NaN',1,'Inf',1)

		if BalancethicknessAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','smb.accumulation','size',[md.mesh.numberofvertices],'NaN',1,'Inf',1)

		if MasstransportAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','smb.runoff','timeseries',1,'NaN',1,'Inf',1)

		if BalancethicknessAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','smb.runoff','size',[md.mesh.numberofvertices],'NaN',1,'Inf',1)

		if MasstransportAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','smb.evaporation','timeseries',1,'NaN',1,'Inf',1)

		if BalancethicknessAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','smb.evaporation','size',[md.mesh.numberofvertices],'NaN',1,'Inf',1)
		
		md = checkfield(md,'fieldname','masstransport.requested_outputs','stringrow',1)

		return md
	# }}}
	def marshall(self,md,fid):    # {{{

		yts=365.0*24.0*3600.0

		WriteData(fid,'enum',SmbEnum(),'data',SMBcomponentsEnum(),'format','Integer');
		WriteData(fid,'object',self,'class','smb','fieldname','accumulation','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'class','smb','fieldname','runoff','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'class','smb','fieldname','evaporation','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		
		#process requested outputs
		outputs = self.requested_outputs
		indices = [i for i, x in enumerate(outputs) if x == 'default']
		if len(indices) > 0:
			outputscopy=outputs[0:max(0,indices[0]-1)]+self.defaultoutputs(md)+outputs[indices[0]+1:]
			outputs    =outputscopy
		WriteData(fid,'data',outputs,'enum',SmbRequestedOutputsEnum(),'format','StringArray')

	# }}}
