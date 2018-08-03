from fielddisplay import fielddisplay
from EnumDefinitions import *
from pairoptions import pairoptions
from checkfield import checkfield
from WriteData import WriteData
from MeshProfileIntersection import MeshProfileIntersection
import os

class massfluxatgate(object):
	"""
	MASSFLUXATEGATE class definition

	   Usage:
		  massfluxatgate=massfluxatgate('GateName','PathToExpFile')
	"""

	def __init__(self,**kwargs): # {{{

		self.name            = ''
		self.definitionenum  = 0
		self.profilename     = ''
		self.segments        = float('NaN')

		#set defaults
		self.setdefaultparameters()

		#use provided options to change fields
		options=pairoptions(**kwargs)

		#OK get other fields
		self=options.AssignObjectFields(self)

		#}}}
	def __repr__(self): # {{{

		string="   Massfluxatgate:"
		string="%s\n%s"%(string,fielddisplay(self,'name','identifier for this massfluxatgate response'))
		string="%s\n%s"%(string,fielddisplay(self,'definitionenum','enum that identifies this output definition uniquely, from Outputdefinition[1-10]Enum'))
		string="%s\n%s"%(string,fielddisplay(self,'profilename','name of file (shapefile or argus file) defining a profile (or gate)'))
		return string
		#}}}
	def setdefaultparameters(self): # {{{
		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{
		
		if  not isinstance(self.name, str):
			raise RuntimeError("massfluxatgate error message: 'name' field should be a string!")
			
		if  not isinstance(self.profilename, str):
			raise RuntimeError("massfluxatgate error message: 'profilename' field should be a string!") 

			md = checkfield(md,'field',self.definitionenum,'values',[Outputdefinition1Enum(),Outputdefinition2Enum(),Outputdefinition3Enum(),Outputdefinition4Enum(),Outputdefinition5Enum(),Outputdefinition6Enum(),Outputdefinition7Enum(),Outputdefinition8Enum(),Outputdefinition9Enum(),Outputdefinition10Enum()])
		
		#check the profilename points to a file!: 
		if not os.path.isfile(self.profilename):
			raise RuntimeError("massfluxatgate error message: file name for profile corresponding to gate does not point to a legitimate file on disk!")

		return md
	# }}}
	def marshall(self,md,fid):    # {{{
		
		#before marshalling, we need to create the segments out of the profilename: 
		self.segments=MeshProfileIntersection(md.mesh.elements,md.mesh.x,md.mesh.y,self.profilename)[0]

		#ok, marshall name and segments: 
		WriteData(fid,'object',self,'fieldname','name','format','String')
		WriteData(fid,'object',self,'fieldname','definitionenum','format','Integer')
		WriteData(fid,'object',self,'fieldname','segments','format','DoubleMat','mattype',1)

	# }}}
