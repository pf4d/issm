import numpy
from project3d import project3d
from collections import OrderedDict
from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData
import MatlabFuncs as m

class qmu(object):
	"""
	QMU class definition

	   Usage:
	      qmu=qmu();
	"""

	def __init__(self): # {{{
		self.isdakota                    = 0
		self.variables                   = OrderedDict()
		self.responses                   = OrderedDict()
		self.method                      = OrderedDict()
		self.params                      = OrderedDict()
		self.results                     = OrderedDict()
		self.partition                   = float('NaN')
		self.numberofpartitions          = 0
		self.numberofresponses           = 0
		self.variabledescriptors         = []
		self.responsedescriptors         = []
		self.mass_flux_profile_directory = float('NaN')
		self.mass_flux_profiles          = float('NaN')
		self.mass_flux_segments          = []
		self.adjacency                   = float('NaN')
		self.vertex_weight               = float('NaN')

		#set defaults
		self.setdefaultparameters()

		#}}}
	def __repr__(self):    # {{{
		s ='   qmu parameters:\n'

		s+="%s\n" % fielddisplay(self,'isdakota','is qmu analysis activated?')
		for i,variable in enumerate(self.variables.items()):
			s+="         variables%s:  (arrays of each variable class)\n" % \
					string_dim(self.variables,i)
			fnames=vars(variable)
			maxlen=0
			for fname in fnames:
				maxlen=max(maxlen,len(fname))

			for fname in fnames:
				s+="'            %-*s:    [%ix%i]    '%s'\n" % \
						(maxlen+1,fname,size(getattr(variable,fname)),type(getattr(variable,fname)))

		for i,response in enumerate(self.responses.items()):
			s+="         responses%s:  (arrays of each response class)\n" % \
					string_dim(self.responses,i)
			fnames=vars(response)
			maxlen=0
			for fname in fnames:
				maxlen=max(maxlen,len(fname))

			for fname in fnames:
				s+="            %-*s:    [%ix%i]    '%s'\n" % \
						(maxlen+1,fname,size(getattr(response,fname)),type(getattr(response,fname)))

		s+="%s\n" % fielddisplay(self,'numberofresponses','number of responses') 

		for i,method in enumerate(self.method.items()):
			if isinstance(method,'dakota_method'):
				s+="            method%s :    '%s'\n" % \
						(string_dim(method,i),method.method)

		for i,param in enumerate(self.params.items()):
			s+="         params%s:  (array of method-independent parameters)\n" % \
					string_dim(self.params,i)
			fnames=vars(param)
			maxlen=0
			for fname in fnames:
				maxlen=max(maxlen,len(fname))

			for fname in fnames:
				s+="            %-*s: %s\n" % \
						(maxlen+1,fname,any2str(getattr(param,fname)))

		for i,result in enumerate(self.results.items()):
			s+="         results%s:  (information from dakota files)\n" % \
					string_dim(self.results,i)
			fnames=vars(result)
			maxlen=0
			for fname in fnames:
				maxlen=max(maxlen,len(fname))

			for fname in fnames:
				s+="            %-*s:    [%ix%i]    '%s'\n" % \
						(maxlen+1,fname,size(getattr(result,fname)),type(getattr(result,fname)))

		s+="%s\n" % fielddisplay(self,'partition','user provided mesh partitioning, defaults to metis if not specified') 
		s+="%s\n" % fielddisplay(self,'numberofpartitions','number of partitions for semi-discrete qmu') 
		s+="%s\n" % fielddisplay(self,'variabledescriptors','')
		s+="%s\n" % fielddisplay(self,'responsedescriptors','')
		s+="%s\n" % fielddisplay(self,'method','array of dakota_method class')
		s+="%s\n" % fielddisplay(self,'mass_flux_profile_directory','directory for mass flux profiles')
		s+="%s\n" % fielddisplay(self,'mass_flux_profiles','list of mass_flux profiles')
		s+="%s\n" % fielddisplay(self,'mass_flux_segments','')
		s+="%s\n" % fielddisplay(self,'adjacency','')
		s+="%s\n" % fielddisplay(self,'vertex_weight','weight applied to each mesh vertex')

		return s
	# }}}
	def extrude(self,md): # {{{
		self.partition=project3d(md,'vector',numpy.transpose(self.partition),'type','node')
		return self
	#}}}
	def setdefaultparameters(self): # {{{
		return self
	#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		#Early return
		if not md.qmu.isdakota:
			return

		if not md.qmu.params.evaluation_concurrency==1:
			md.checkmessage("concurrency should be set to 1 when running dakota in library mode")
		if md.qmu.partition:
			if not numpy.size(md.qmu.partition)==md.mesh.numberofvertices:
				md.checkmessage("user supplied partition for qmu analysis should have size md.mesh.numberofvertices x 1")
			if not min(md.qmu.partition)==0:
				md.checkmessage("partition vector not indexed from 0 on")
			if max(md.qmu.partition)>=md.qmu.numberofpartitions:
				md.checkmessage("for qmu analysis, partitioning vector cannot go over npart, number of partition areas")

		if not m.strcmpi(md.cluster.name,'none'):
			if not md.settings.waitonlock:
				md.checkmessage("waitonlock should be activated when running qmu in parallel mode!")

		return md
	# }}}
	def marshall(self,md,fid):    # {{{
		WriteData(fid,'object',self,'fieldname','isdakota','format','Boolean')
		if not self.isdakota:
			WriteData(fid,'data',False,'enum',QmuMassFluxSegmentsPresentEnum(),'format','Boolean');
			return
		WriteData(fid,'object',self,'fieldname','partition','format','DoubleMat','mattype',2)
		WriteData(fid,'object',self,'fieldname','numberofpartitions','format','Integer')
		WriteData(fid,'object',self,'fieldname','numberofresponses','format','Integer')
		WriteData(fid,'object',self,'fieldname','variabledescriptors','format','StringArray')
		WriteData(fid,'object',self,'fieldname','responsedescriptors','format','StringArray')
		if not self.mass_flux_segments:
			WriteData(fid,'data',self.mass_flux_segments,'enum',MassFluxSegmentsEnum(),'format','MatArray');
			flag=True; 
		else:
			flag=False; 
		WriteData(fid,'data',flag,'enum',QmuMassFluxSegmentsPresentEnum(),'format','Boolean');
	# }}}
