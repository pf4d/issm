import numpy
from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData
from project3d import project3d

class SMBforcing(object):
	"""
	SMBforcing Class definition

	   Usage:
	      SMB=SMBforcing();
	"""

	def __init__(self): # {{{
		self.mass_balance = float('NaN')
		self.requested_outputs      = []
		#}}}
	def __repr__(self): # {{{
		string="   surface forcings parameters:"
		string="%s\n%s"%(string,fielddisplay(self,'mass_balance','surface mass balance [m/yr ice eq]'))
		string="%s\n%s"%(string,fielddisplay(self,'requested_outputs','additional outputs requested'))
		return string
		#}}}
	def extrude(self,md): # {{{

		self.mass_balance=project3d(md,'vector',self.mass_balance,'type','node');
		return self
	#}}}
	def defaultoutputs(self,md): # {{{
		return []
	#}}}
	def initialize(self,md): # {{{

		if numpy.all(numpy.isnan(self.mass_balance)):
			self.mass_balance=numpy.zeros((md.mesh.numberofvertices,1))
			print("      no SMBforcing.mass_balance specified: values set as zero")

		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		if MasstransportAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','smb.mass_balance','timeseries',1,'NaN',1,'Inf',1)

		if BalancethicknessAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','smb.mass_balance','size',[md.mesh.numberofvertices],'NaN',1,'Inf',1)

		md = checkfield(md,'fieldname','masstransport.requested_outputs','stringrow',1)
		return md
	# }}}
	def marshall(self,md,fid):    # {{{

		yts=365.0*24.0*3600.0

		WriteData(fid,'enum',SmbEnum(),'data',SMBforcingEnum(),'format','Integer');
		WriteData(fid,'object',self,'class','smb','fieldname','mass_balance','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		
		#process requested outputs
		outputs = self.requested_outputs
		indices = [i for i, x in enumerate(outputs) if x == 'default']
		if len(indices) > 0:
			outputscopy=outputs[0:max(0,indices[0]-1)]+self.defaultoutputs(md)+outputs[indices[0]+1:]
			outputs    =outputscopy
		WriteData(fid,'data',outputs,'enum',SmbRequestedOutputsEnum(),'format','StringArray')

	# }}}
