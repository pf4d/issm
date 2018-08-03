from fielddisplay import fielddisplay
from EnumDefinitions import *
from StringToEnum import StringToEnum
from checkfield import checkfield
from WriteData import WriteData
import numpy as np

class outputdefinition(object):
	"""
	OUTPUTDEFINITION class definition

	   Usage:
	      outputdefinition=outputdefinition();
	"""

	def __init__(self): # {{{
		self.definitions                   = []
		#}}}
	def __repr__(self): # {{{
		string="   Outputdefinitions:"

		string="%s\n%s"%(string,fielddisplay(self,"definitions","list of potential outputs that can be requested, but which need additional data to be defined"))

		return string
		#}}}
	def setdefaultparameters(self): # {{{
		return self
		#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{
		
		md = checkfield(md,'fieldname','outputdefinition.definitions','cell',1)
		for definition in self.definitions:
			definition.checkconsistency(md,solution,analyses);

	# }}}
	def marshall(self,md,fid):    # {{{
		
		enums=np.zeros(len(self.definitions),)
		
		for i in range(len(self.definitions)):
			self.definitions[i].marshall(md,fid);
			classdefinition=self.definitions[i].__class__.__name__
			classdefinition=classdefinition[0].upper()+classdefinition[1:]
			enums[i]=StringToEnum(classdefinition)[0]
		
		enums=np.unique(enums);
		
		WriteData(fid,'data',enums,'enum',OutputdefinitionListEnum(),'format','DoubleMat','mattype',1);
	# }}}
