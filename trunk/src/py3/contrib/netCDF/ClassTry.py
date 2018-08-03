#module imports {{{
from netCDF4 import Dataset
import time
import collections
from os import path, remove
#}}}

				
class truc(object):
	#properties
	def __init__(self,*filename):#{{{

		def netCDFread(filename):
			def walktree(data):
				keys = list(data.groups.keys())
				yield keys
				for key in keys:
					for children in walktree(data.groups[str(key)]):
						yield children

			if path.exists(filename):
				print(('Opening {} for reading '.format(filename)))
				NCData=Dataset(filename, 'r')
				class_dict={}
				
				for children in walktree(NCData):
					for child in children:
						class_dict[str(child)]=str(getattr(NCData.groups[str(child)],'classtype'))

				return class_dict

		if filename:		
			classtype=netCDFread(filename[0])
		else:
			classtype=self.default_prop()
			
		module=list(map(__import__,dict.values(classtype)))

		for i,mod in enumerate(dict.keys(classtype)):
			self.__dict__[mod] = getattr(module[i],str(classtype[str(mod)]))()
			
		#}}}
	def default_prop(self):    # {{{
		# ordered list of properties since vars(self) is random
		return {'mesh':'mesh2d',\
		        'mask':'mask',\
		        'geometry':'geometry',\
		        'constants':'constants',\
		        'smb':'SMB',\
		        'basalforcings':'basalforcings',\
		        'materials':'matice',\
		        'damage':'damage',\
		        'friction':'friction',\
		        'flowequation':'flowequation',\
		        'timestepping':'timestepping',\
		        'initialization':'initialization',\
		        'rifts':'rifts',\
		        'debug':'debug',\
		        'verbose':'verbose',\
		        'settings':'settings',\
		        'toolkits':'toolkits',\
		        'cluster':'generic',\
		        'balancethickness':'balancethickness',\
		        'stressbalance':'stressbalance',\
		        'groundingline':'groundingline',\
		        'hydrology':'hydrologyshreve',\
		        'masstransport':'masstransport',\
		        'thermal':'thermal',\
		        'steadystate':'steadystate',\
		        'transient':'transient',\
		        'calving':'calving',\
						'gia':'gia',\
		        'autodiff':'autodiff',\
		        'flaim':'flaim',\
		        'inversion':'inversion',\
		        'qmu':'qmu',\
		        'outputdefinition':'outputdefinition',\
		        'results':'results',\
		        'radaroverlay':'radaroverlay',\
		        'miscellaneous':'miscellaneous',\
		        'private':'private'}
	# }}}
		
	def __repr__(obj): #{{{
		#print "Here %s the number: %d" % ("is", 37)
		string="%19s: %-22s -- %s" % ("mesh","[%s,%s]" % ("1x1",obj.mesh.__class__.__name__),"mesh properties")
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("mask","[%s,%s]" % ("1x1",obj.mask.__class__.__name__),"defines grounded and floating elements"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("geometry","[%s,%s]" % ("1x1",obj.geometry.__class__.__name__),"surface elevation, bedrock topography, ice thickness,..."))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("constants","[%s,%s]" % ("1x1",obj.constants.__class__.__name__),"physical constants"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("smb","[%s,%s]" % ("1x1",obj.smb.__class__.__name__),"surface forcings"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("basalforcings","[%s,%s]" % ("1x1",obj.basalforcings.__class__.__name__),"bed forcings"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("materials","[%s,%s]" % ("1x1",obj.materials.__class__.__name__),"material properties"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("damage","[%s,%s]" % ("1x1",obj.damage.__class__.__name__),"damage propagation laws"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("friction","[%s,%s]" % ("1x1",obj.friction.__class__.__name__),"basal friction/drag properties"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("flowequation","[%s,%s]" % ("1x1",obj.flowequation.__class__.__name__),"flow equations"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("timestepping","[%s,%s]" % ("1x1",obj.timestepping.__class__.__name__),"time stepping for transient models"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("initialization","[%s,%s]" % ("1x1",obj.initialization.__class__.__name__),"initial guess/state"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("rifts","[%s,%s]" % ("1x1",obj.rifts.__class__.__name__),"rifts properties"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("debug","[%s,%s]" % ("1x1",obj.debug.__class__.__name__),"debugging tools (valgrind, gprof)"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("verbose","[%s,%s]" % ("1x1",obj.verbose.__class__.__name__),"verbosity level in solve"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("settings","[%s,%s]" % ("1x1",obj.settings.__class__.__name__),"settings properties"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("toolkits","[%s,%s]" % ("1x1",obj.toolkits.__class__.__name__),"PETSc options for each solution"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("cluster","[%s,%s]" % ("1x1",obj.cluster.__class__.__name__),"cluster parameters (number of cpus...)"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("balancethickness","[%s,%s]" % ("1x1",obj.balancethickness.__class__.__name__),"parameters for balancethickness solution"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("stressbalance","[%s,%s]" % ("1x1",obj.stressbalance.__class__.__name__),"parameters for stressbalance solution"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("groundingline","[%s,%s]" % ("1x1",obj.groundingline.__class__.__name__),"parameters for groundingline solution"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("hydrology","[%s,%s]" % ("1x1",obj.hydrology.__class__.__name__),"parameters for hydrology solution"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("masstransport","[%s,%s]" % ("1x1",obj.masstransport.__class__.__name__),"parameters for masstransport solution"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("thermal","[%s,%s]" % ("1x1",obj.thermal.__class__.__name__),"parameters for thermal solution"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("steadystate","[%s,%s]" % ("1x1",obj.steadystate.__class__.__name__),"parameters for steadystate solution"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("transient","[%s,%s]" % ("1x1",obj.transient.__class__.__name__),"parameters for transient solution"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("calving","[%s,%s]" % ("1x1",obj.calving.__class__.__name__),"parameters for calving"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("autodiff","[%s,%s]" % ("1x1",obj.autodiff.__class__.__name__),"automatic differentiation parameters"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("flaim","[%s,%s]" % ("1x1",obj.flaim.__class__.__name__),"flaim parameters"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("inversion","[%s,%s]" % ("1x1",obj.inversion.__class__.__name__),"parameters for inverse methods"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("qmu","[%s,%s]" % ("1x1",obj.qmu.__class__.__name__),"dakota properties"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("outputdefinition","[%s,%s]" % ("1x1",obj.outputdefinition.__class__.__name__),"output definition"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("results","[%s,%s]" % ("1x1",obj.results.__class__.__name__),"model results"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("radaroverlay","[%s,%s]" % ("1x1",obj.radaroverlay.__class__.__name__),"radar image for plot overlay"))
		string="%s\n%s" % (string,"%19s: %-22s -- %s" % ("miscellaneous","[%s,%s]" % ("1x1",obj.miscellaneous.__class__.__name__),"miscellaneous fields"))
		return string
	# }}}
