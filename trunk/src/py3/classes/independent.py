import numpy
from pairoptions import pairoptions
from fielddisplay import fielddisplay
import MatlabFuncs as m
from EnumDefinitions import *

class independent(object):
	"""
	INDEPENDENT class definition

	   Usage:
	      independent=independent();
	"""

	def __init__(self,**kwargs):    # {{{
		self.name                 = ''
		self.type                 = ''
		self.fos_forward_index    = float('NaN')
		self.fov_forward_indices  = numpy.array([])
		self.nods                 = 0

		#set defaults
		self.setdefaultparameters()

		#use provided options to change fields
		options=pairoptions(**kwargs)

		#OK get other fields
		self=options.AssignObjectFields(self)
	# }}}
	def __repr__(self):    # {{{
		s ="   independent variable:\n"

		s+="%s\n" % fielddisplay(self,'name',"variable name (must match corresponding Enum)")
		s+="%s\n" % fielddisplay(self,'type',"type of variable ('vertex' or 'scalar')")
		if not numpy.isnan(self.fos_forward_index):
			s+="%s\n" % fielddisplay(self,'fos_forward_index',"index for fos_foward driver of ADOLC")
		if numpy.any(numpy.logical_not(numpy.isnan(self.fov_forward_indices))):
			s+="%s\n" % fielddisplay(self,'fov_forward_indices',"indices for fov_foward driver of ADOLC")

		return s
	# }}}
	def setdefaultparameters(self):    # {{{
		#do nothing
		return self
	# }}}
	def checkconsistency(self,md,i,solution,analyses,driver):    # {{{
		if not numpy.isnan(self.fos_forward_index):
			if not strcmpi(driver,'fos_forward'):
				raise TypeError("cannot declare an independent with a fos_forward_index when the driver is not fos_forward!")
			if self.nods==0:
				raise TypeError("independent checkconsistency error: nods should be set to the size of the independent variable")

		if self.fov_forward_indices:
			if not strcmpi(driver,'fov_forward'):
				raise TypeError("cannot declare an independent with fov_forward_indices when the driver is not fov_forward!")
			if self.nods==0:
				raise TypeError("independent checkconsistency error: nods should be set to the size of the independent variable")
			md = checkfield(md,'fieldname',"autodiff.independents[%d].fov_forward_indices" % i,'>=',1,'<=',self.nods,'size',[float('NaN'),1])

		return md
	# }}}
	def typetoscalar(self):    # {{{
		if   strcmpi(self.type,'scalar'):
			scalar=0
		elif strcmpi(self.type,'vertex'):
			scalar=1

		return scalar
	# }}}
