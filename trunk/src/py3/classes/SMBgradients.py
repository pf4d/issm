from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData
from project3d import project3d

class SMBgradients(object):
	"""
	SMBgradients Class definition

	   Usage:
	      SMBgradients=SMBgradients();
	"""

	def __init__(self): # {{{
		self.href    = float('NaN')
		self.smbref  = float('NaN')
		self.b_pos   = float('NaN')
		self.b_neg   = float('NaN')
		self.requested_outputs      = []
		#}}}
	def __repr__(self): # {{{
		string="   surface forcings parameters:"

		string="%s\n%s"%(string,fielddisplay(self,'issmbgradients','is smb gradients method activated (0 or 1, default is 0)'))
		string="%s\n%s"%(string,fielddisplay(self,'href',' reference elevation from which deviation is used to calculate SMB adjustment in smb gradients method'))
		string="%s\n%s"%(string,fielddisplay(self,'smbref',' reference smb from which deviation is calculated in smb gradients method'))
		string="%s\n%s"%(string,fielddisplay(self,'b_pos',' slope of hs - smb regression line for accumulation regime required if smb gradients is activated'))
		string="%s\n%s"%(string,fielddisplay(self,'b_neg',' slope of hs - smb regression line for ablation regime required if smb gradients is activated'))
		string="%s\n%s"%(string,fielddisplay(self,'requested_outputs','additional outputs requested'))

		return string
		#}}}
	def extrude(self,md): # {{{

		#Nothing for now
		return self
	#}}}
	def defaultoutputs(self,md): # {{{
		return []
	#}}}
	def initialize(self,md): # {{{

		#Nothing for now

		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		if MasstransportAnalysisEnum() in analyses:
			md = checkfield(md,'fieldname','smb.href','timeseries',1,'NaN',1,'Inf',1)
			md = checkfield(md,'fieldname','smb.smbref','timeseries',1,'NaN',1,'Inf',1)
			md = checkfield(md,'fieldname','smb.b_pos','timeseries',1,'NaN',1,'Inf',1)
			md = checkfield(md,'fieldname','smb.b_neg','timeseries',1,'NaN',1,'Inf',1)

		md = checkfield(md,'fieldname','masstransport.requested_outputs','stringrow',1)
		return md
	# }}}
	def marshall(self,md,fid):    # {{{

		yts=365.0*24.0*3600.0

		WriteData(fid,'enum',SmbEnum(),'data',SMBgradientsEnum(),'format','Integer');
		WriteData(fid,'object',self,'class','smb','fieldname','href','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'class','smb','fieldname','smbref','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'class','smb','fieldname','b_pos','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		WriteData(fid,'object',self,'class','smb','fieldname','b_neg','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
		
		#process requested outputs
		outputs = self.requested_outputs
		indices = [i for i, x in enumerate(outputs) if x == 'default']
		if len(indices) > 0:
			outputscopy=outputs[0:max(0,indices[0]-1)]+self.defaultoutputs(md)+outputs[indices[0]+1:]
			outputs    =outputscopy
		WriteData(fid,'data',outputs,'enum',SmbRequestedOutputsEnum(),'format','StringArray')

	# }}}
