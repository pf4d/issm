from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData

class transient(object):
	"""
	TRANSIENT class definition

	   Usage:
	      transient=transient();
	"""

	def __init__(self): # {{{
		self.issmb   = False
		self.ismasstransport   = False
		self.isstressbalance   = False
		self.isthermal         = False
		self.isgroundingline   = False
		self.isgia             = False
		self.isdamageevolution = False
		self.islevelset        = False
		self.iscalving         = False
		self.ishydrology       = False
		self.requested_outputs = []

		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self): # {{{
		string='   transient solution parameters:'
		string="%s\n%s"%(string,fielddisplay(self,'issmb','indicates if a surface mass balance solution is used in the transient'))
		string="%s\n%s"%(string,fielddisplay(self,'ismasstransport','indicates if a masstransport solution is used in the transient'))
		string="%s\n%s"%(string,fielddisplay(self,'isstressbalance','indicates if a stressbalance solution is used in the transient'))
		string="%s\n%s"%(string,fielddisplay(self,'isthermal','indicates if a thermal solution is used in the transient'))
		string="%s\n%s"%(string,fielddisplay(self,'isgroundingline','indicates if a groundingline migration is used in the transient'))
		string="%s\n%s"%(string,fielddisplay(self,'isgia','indicates if a postglacial rebound is used in the transient'))
		string="%s\n%s"%(string,fielddisplay(self,'isdamageevolution','indicates whether damage evolution is used in the transient'))
		string="%s\n%s"%(string,fielddisplay(self,'islevelset','LEVELSET METHOD DESCRIPTION'))
		string="%s\n%s"%(string,fielddisplay(self,'iscalving','indicates whether calving is used in the transient'))
		string="%s\n%s"%(string,fielddisplay(self,'ishydrology','indicates whether an hydrology model is used'))
		string="%s\n%s"%(string,fielddisplay(self,'requested_outputs','list of additional outputs requested'))
		return string
		#}}}
	def defaultoutputs(self,md): # {{{

		if self.issmb:
			return ['SmbMassBalance']
		else:
			return []

	#}}}
	def setallnullparameters(self): # {{{
		
		#Nothing done
		self.issmb   = False
		self.ismasstransport   = False
		self.isstressbalance   = False
		self.isthermal         = False
		self.isgroundingline   = False
		self.isgia             = False
		self.isdamageevolution = False
		self.islevelset        = False
		self.iscalving         = False
		self.ishydrology       = False

		#default output
		self.requested_outputs=[]
		return self
	#}}}
	def setdefaultparameters(self): # {{{
		
		#full analysis: Stressbalance, Masstransport and Thermal but no groundingline migration for now
		self.issmb = True
		self.ismasstransport = True
		self.isstressbalance = True
		self.isthermal       = True
		self.isgroundingline = False
		self.isgia           = False
		self.isdamageevolution = False
		self.islevelset      = False
		self.iscalving       = False
		self.ishydrology     = False

		#default output
		self.requested_outputs=['default']
		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		#Early return
		if not solution==TransientSolutionEnum():
			return md

		md = checkfield(md,'fieldname','transient.issmb','numel',[1],'values',[0,1])
		md = checkfield(md,'fieldname','transient.ismasstransport','numel',[1],'values',[0,1])
		md = checkfield(md,'fieldname','transient.isstressbalance','numel',[1],'values',[0,1])
		md = checkfield(md,'fieldname','transient.isthermal','numel',[1],'values',[0,1])
		md = checkfield(md,'fieldname','transient.isgroundingline','numel',[1],'values',[0,1])
		md = checkfield(md,'fieldname','transient.isgia','numel',[1],'values',[0,1])
		md = checkfield(md,'fieldname','transient.isdamageevolution','numel',[1],'values',[0,1])
		md = checkfield(md,'fieldname','transient.islevelset','numel',[1],'values',[0,1])
		md = checkfield(md,'fieldname','transient.ishydrology','numel',[1],'values',[0,1])
		md = checkfield(md,'fieldname','transient.iscalving','numel',[1],'values',[0,1]);
		md = checkfield(md,'fieldname','transient.requested_outputs','stringrow',1)

		return md
	# }}}
	def marshall(self,md,fid):    # {{{
		WriteData(fid,'object',self,'fieldname','issmb','format','Boolean')
		WriteData(fid,'object',self,'fieldname','ismasstransport','format','Boolean')
		WriteData(fid,'object',self,'fieldname','isstressbalance','format','Boolean')
		WriteData(fid,'object',self,'fieldname','isthermal','format','Boolean')
		WriteData(fid,'object',self,'fieldname','isgroundingline','format','Boolean')
		WriteData(fid,'object',self,'fieldname','isgia','format','Boolean')
		WriteData(fid,'object',self,'fieldname','isdamageevolution','format','Boolean')
		WriteData(fid,'object',self,'fieldname','islevelset','format','Boolean')
		WriteData(fid,'object',self,'fieldname','ishydrology','format','Boolean')
		WriteData(fid,'object',self,'fieldname','iscalving','format','Boolean')

		#process requested outputs
		outputs = self.requested_outputs
		indices = [i for i, x in enumerate(outputs) if x == 'default']
		if len(indices) > 0:
			outputscopy=outputs[0:max(0,indices[0]-1)]+self.defaultoutputs(md)+outputs[indices[0]+1:]
			outputs    =outputscopy
		WriteData(fid,'data',outputs,'enum',TransientRequestedOutputsEnum(),'format','StringArray')
	# }}}
