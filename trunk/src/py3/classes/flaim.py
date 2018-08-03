import numpy
from collections import OrderedDict
from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData

class flaim(object):
	"""
	FLAIM class definition

	   Usage:
	      flaim=flaim();
	"""

	def __init__(self): # {{{
		self.targets            = ''
		self.tracks             = ''
		self.flightreqs         = OrderedDict()
		self.criterion          = float('NaN')
		self.gridsatequator     = 200000
		self.usevalueordering   = True
		self.split_antimeridian = True
		self.solution           = ''
		self.quality            = 0
		self.path_optimize      = False
		self.opt_ndir           = 1
		self.opt_dist           = 25
		self.opt_niter          = 30000
		#}}}
	def __repr__(self): # {{{
		string='   FLAIM - Flight Line Adaptation using Ice sheet Modeling:'

		string="%s\n\n%s"%(string,'      Input:')
		string="%s\n%s"%(string,fielddisplay(self,'targets'            ,'name of kml output targets file '))
		string="%s\n%s"%(string,fielddisplay(self,'tracks'             ,'name of kml input tracks file '))
		string="%s\n%s"%(string,fielddisplay(self,'flightreqs'         ,'structure of kml flight requirements (not used yet)'))
		string="%s\n%s"%(string,fielddisplay(self,'criterion'          ,'element or nodal criterion for flight path evaluation (metric)'))

		string="%s\n\n%s"%(string,'      Arguments:')
		string="%s\n%s"%(string,fielddisplay(self,'gridsatequator'     ,'number of grids at equator (determines resolution)'))
		string="%s\n%s"%(string,fielddisplay(self,'usevalueordering'   ,'flag to consider target values for flight path evaluation'))
		string="%s\n%s"%(string,fielddisplay(self,'split_antimeridian' ,'flag to split polygons on the antimeridian'))
		
		string="%s\n\n%s"%(string,'      Optimization:')
		string="%s\n%s"%(string,fielddisplay(self,'path_optimize'     ,'optimize? (default false)'))
		string="%s\n%s"%(string,fielddisplay(self,'opt_ndir'     ,['number of directions to test when moving a point.  If this value = 1, a random direction is tested.',\
										  'A value > 1 results in directions equally spaced from [0, 2*PI] being tested.',\
										  'For example, 4 would result in directions [0, PI/2, PI, 3PI/2].']))
		string="%s\n%s"%(string,fielddisplay(self,'opt_dist'     ,'specifies the distance in km (default 25) to move a randomly selected path point on each iteration'))
		string="%s\n%s"%(string,fielddisplay(self,'opt_niter'     ,['number of iterations (default 30,000) to run for flightplan optimization',\
										   'i.e. the number of times to randomly select a point and move it.']))

		string="%s\n\n%s"%(string,'      Output:')
		string="%s\n%s"%(string,fielddisplay(self,'solution'           ,'name of kml solution file'))
		string="%s\n%s"%(string,fielddisplay(self,'quality'            ,'quality of kml solution'))
		return string
		#}}}
	def checkconsistency(self,md,solution,analyses):    # {{{

		#Early return
		if not solution==FlaimSolutionEnum():
			return md

		md = checkfield(md,'fieldname','flaim.tracks','file',1)
		if numpy.any(numpy.isnan(md.flaim.criterion)) or not md.flaim.criterion:
			md = checkfield(md,'fieldname','flaim.targets','file',1)
		else:
			md = checkfield(md,'fieldname','flaim.criterion','numel',[md.mesh.numberofvertices,md.mesh.numberofelements])

		return md
	# }}}
