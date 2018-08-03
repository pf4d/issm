"""
== == == == == == == == == == == == == == == == == == ==
Auto generated python script for ISSM:   /home/andrei/issm/trunk-jpl/src/m/classes/adinversion.m
Created on 2015-05-15 via translateToPy.py Ver 1.0 by andrei
== == == == == == == == == == == == == == == == == == ==

Matlab script conversion into python
translateToPy.py Author: Michael Pellegrin
translateToPy.py Date: 09/24/12
== == == == == == == == == == == == == == == == == == ==
"""

from MatlabFuncs import *

from EnumDefinitions import *
from numpy import *

# ADINVERSION class definition

# 

#    Usage:

#       adinversion=adinversion();



class adinversion:
	def __init__(self): 
		iscontrol                   = 0
		control_parameters          = float('Nan')
		control_scaling_factors     = float('Nan')
		maxsteps                    = 0
		maxiter                     = 0
		dxmin                       = 0
		gttol                       = 0
		cost_functions              = float('Nan')
		cost_functions_coefficients = float('Nan')
		min_parameters              = float('Nan')
		max_parameters              = float('Nan')
		vx_obs                      = float('Nan')
		vy_obs                      = float('Nan')
		vz_obs                      = float('Nan')
		vel_obs                     = float('Nan')
		thickness_obs               = float('Nan')
		surface_obs                 = float('Nan')

	def setdefaultparameters(self):

		self.control_parameters=['FrictionCoefficient']


# 		Scaling factor for each control
		self.control_scaling_factors=1

# 		number of iterations
		self.maxsteps=20
		self.maxiter=40

#		several responses can be used:
		self.cost_functions=['FrictionCoefficient']

# 		m1qn3 parameters
		self.dxmin  = 0.1
		self.gttol = 1e-4

		return self
	
	def checkconsistency(self, md, solution, analyses): 

# 			Early return
		if not self.iscontrol:
			return

		if not IssmConfig('_HAVE_M1QN3_'):
			md = checkmessage(md,['M1QN3 has not been installed, ISSM needs to be reconfigured and recompiled with AD'])


		num_controls=numpy.numel(md.inversion.control_parameters)
		num_costfunc=numpy.size(md.inversion.cost_functions,2)


		md = checkfield(md,'fieldname','inversion.iscontrol','values',[0, 1])
		md = checkfield(md,'fieldname','inversion.control_parameters','cell',1,'values',\
			['BalancethicknessThickeningRate' 'FrictionCoefficient' 'MaterialsRheologyBbar' 'DamageDbar',\
			'Vx' 'Vy' 'Thickness' 'BalancethicknessOmega' 'BalancethicknessApparentMassbalance'])
		md = checkfield(md,'fieldname','inversion.control_scaling_factors','size',[1, num_controls],'>',0,float('Nan'),1)
		md = checkfield(md,'fieldname','inversion.maxsteps','numel',1,'>=',0)
		md = checkfield(md,'fieldname','inversion.maxiter','numel',1,'>=',0)
		md = checkfield(md,'fieldname','inversion.dxmin','numel',1,'>',0)
		md = checkfield(md,'fieldname','inversion.gttol','numel',1,'>',0)
		md = checkfield(md,'fieldname','inversion.cost_functions','size',[1, num_costfunc],'values', [i for i in range(101,106)]+[201]+[i for i in range(501,507)]+[i for i in range(601,605)]+[i for i in range(1001, 1011)])
		md = checkfield(md,'fieldname','inversion.cost_functions_coefficients','size',[md.mesh.numberofvertices, num_costfunc],'>=',0)
		md = checkfield(md,'fieldname','inversion.min_parameters','size',[md.mesh.numberofvertices, num_controls])
		md = checkfield(md,'fieldname','inversion.max_parameters','size',[md.mesh.numberofvertices, num_controls])


		if solution==BalancethicknessSolutionEnum():
			md = checkfield(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices, 1],float('Nan'),1)
			md = checkfield(md,'fieldname','inversion.surface_obs','size',[md.mesh.numberofvertices, 1], float('Nan'),1)
		elif solution==BalancethicknessSoftSolutionEnum():
			md = checkfield(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices, 1],float('Nan'),1)
		else:
			md = checkfield(md,'fieldname','inversion.vx_obs','size',[md.mesh.numberofvertices, 1],float('Nan'),1)
			if not numpy.strcmp(domaintype(md.mesh),'2Dvertical'):
				md = checkfield(md,'fieldname','inversion.vy_obs','size',[md.mesh.numberofvertices, 1],float('Nan'),1)
		return md

	def __repr__(self):
		string = '   adinversion parameters:'
		string ="%s\n\%s"%(string, fielddisplay(self,'iscontrol','is inversion activated?'))
		string ="%s\n\%s"%(string, fielddisplay(self,'control_parameters','ex: [''FrictionCoefficient''], or [''MaterialsRheologyBbar'']'))
		string ="%s\n\%s"%(string, fielddisplay(self,'control_scaling_factors','order of magnitude of each control (useful for multi-parameter optimization)'))
		string ="%s\n\%s"%(string, fielddisplay(self,'maxsteps','maximum number of iterations (gradient computation)'))
		string ="%s\n\%s"%(string, fielddisplay(self,'maxiter','maximum number of Function evaluation (forward run)'))
		string ="%s\n\%s"%(string, fielddisplay(self,'dxmin','convergence criterion: two points less than dxmin from eachother (sup-norm) are considered identical'))
		string ="%s\n\%s"%(string, fielddisplay(self,'gttol','convergence criterion: ||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)'))
		string ="%s\n\%s"%(string, fielddisplay(self,'cost_functions','indicate the type of response for each optimization step'))
		string ="%s\n\%s"%(string, fielddisplay(self,'cost_functions_coefficients','cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter'))
		string ="%s\n\%s"%(string, fielddisplay(self,'min_parameters','absolute minimum acceptable value of the inversed parameter on each vertex'))
		string ="%s\n\%s"%(string, fielddisplay(self,'max_parameters','absolute maximum acceptable value of the inversed parameter on each vertex'))
		string ="%s\n\%s"%(string, fielddisplay(self,'vx_obs','observed velocity x component [m/yr]'))
		string ="%s\n\%s"%(string, fielddisplay(self,'vy_obs','observed velocity y component [m/yr]'))
		string ="%s\n\%s"%(string, fielddisplay(self,'vel_obs','observed velocity magnitude [m/yr]'))
		string ="%s\n\%s"%(string, fielddisplay(self,'thickness_obs','observed thickness [m]'))
		string ="%s\n\%s"%(string, fielddisplay(self,'surface_obs','observed surface elevation [m]'))
		string ="%s\n%s"%(string,'Available cost functions:');
		string ="%s\n%s"%(string,'   101: SurfaceAbsVelMisfit');
		string ="%s\n%s"%(string,'   102: SurfaceRelVelMisfit');
		string ="%s\n%s"%(string,'   103: SurfaceLogVelMisfit');
		string ="%s\n%s"%(string,'   104: SurfaceLogVxVyMisfit');
		string ="%s\n%s"%(string,'   105: SurfaceAverageVelMisfit');
		string ="%s\n%s"%(string,'   201: ThicknessAbsMisfit');
		string ="%s\n%s"%(string,'   501: DragCoefficientAbsGradient');
		string ="%s\n%s"%(string,'   502: RheologyBbarAbsGradient');
		string ="%s\n%s"%(string,'   503: ThicknessAbsGradient');
		
		return string

	def marshall(self):

		yts=365.0*24.0*3600.0;

		WriteData(fid,'object',self,'class','inversion','fieldname','iscontrol','format','Boolean');
		WriteData(fid,'enum',InversionTypeEnum(),'data',4,'format','Integer');
		if not self.iscontrol:
			return
		WriteData(fid,'object',self,'class','inversion','fieldname','control_scaling_factors','format','DoubleMat','mattype',3);
		WriteData(fid,'object',self,'class','inversion','fieldname','maxsteps','format','Integer');
		WriteData(fid,'object',self,'class','inversion','fieldname','maxiter','format','Integer');
		WriteData(fid,'object',self,'class','inversion','fieldname','dxmin','format','Double');
		WriteData(fid,'object',self,'class','inversion','fieldname','gttol','format','Double');
		WriteData(fid,'object',self,'class','inversion','fieldname','cost_functions_coefficients','format','DoubleMat','mattype',1);
		WriteData(fid,'object',self,'class','inversion','fieldname','min_parameters','format','DoubleMat','mattype',3);
		WriteData(fid,'object',self,'class','inversion','fieldname','max_parameters','format','DoubleMat','mattype',3);
		WriteData(fid,'object',self,'class','inversion','fieldname','vx_obs','format','DoubleMat','mattype',1,'scale',1./yts);
		WriteData(fid,'object',self,'class','inversion','fieldname','vy_obs','format','DoubleMat','mattype',1,'scale',1./yts);
		WriteData(fid,'object',self,'class','inversion','fieldname','vz_obs','format','DoubleMat','mattype',1,'scale',1./yts);
		if(numel(self.thickness_obs)==md.mesh.numberofelements):
			mattype=2;
		else:
			mattype=1;
		
		WriteData(fid,'object',self,'class','inversion','fieldname','thickness_obs','format','DoubleMat','mattype',mattype);
		WriteData(fid,'object',self,'class','inversion','fieldname','surface_obs','format','DoubleMat','mattype',mattype);

		#process control parameters
		num_control_parameters = numpy.numel(self.control_parameters);
		data = numpy.array([StringToEnum(self.control_parameter[0]) for control_parameter in self.control_parameters]).reshape(1,-1)

		WriteData(fid,'data',data,'enum',InversionControlParametersEnum(),'format','DoubleMat','mattype',3);
		WriteData(fid,'data',num_control_parameters,'enum',InversionNumControlParametersEnum(),'format','Integer');

		#process cost functions
		num_cost_functions=numpy.size(self.cost_functions,2);
		data=copy.deepcopy(self.cost_functions)
		data[numpy.nonzero(self.cost_functions==101)] =SurfaceAbsVelMisfitEnum();
		data[numpy.nonzero(self.cost_functions==102)]=SurfaceRelVelMisfitEnum();
		data[numpy.nonzero(self.cost_functions==103)]=SurfaceLogVelMisfitEnum();
		data[numpy.nonzero(self.cost_functions==104)]=SurfaceLogVxVyMisfitEnum();
		data[numpy.nonzero(self.cost_functions==105)]=SurfaceAverageVelMisfitEnum();
		data[numpy.nonzero(self.cost_functions==201)]=ThicknessAbsMisfitEnum();
		data[numpy.nonzero(self.cost_functions==501)]=DragCoefficientAbsGradientEnum();
		data[numpy.nonzero(self.cost_functions==502)]=RheologyBbarAbsGradientEnum();
		data[numpy.nonzero(self.cost_functions==503)]=ThicknessAbsGradientEnum();
		data[numpy.nonzero(self.cost_functions==504)]=ThicknessAlongGradientEnum();
		data[numpy.nonzero(self.cost_functions==505)]=ThicknessAcrossGradientEnum();
		data[numpy.nonzero(self.cost_functions==506)]=BalancethicknessMisfitEnum();
		data[numpy.nonzero(self.cost_functions==601)]=SurfaceAbsMisfitEnum();
		data[numpy.nonzero(self.cost_functions==1001)]=Outputdefinition1Enum();
		data[numpy.nonzero(self.cost_functions==1002)]=Outputdefinition2Enum();
		data[numpy.nonzero(self.cost_functions==1003)]=Outputdefinition3Enum();
		data[numpy.nonzero(self.cost_functions==1004)]=Outputdefinition4Enum();
		data[numpy.nonzero(self.cost_functions==1005)]=Outputdefinition5Enum();
		data[numpy.nonzero(self.cost_functions==1006)]=Outputdefinition6Enum();
		data[numpy.nonzero(self.cost_functions==1007)]=Outputdefinition7Enum();
		data[numpy.nonzero(self.cost_functions==1008)]=Outputdefinition8Enum();
		data[numpy.nonzero(self.cost_functions==1009)]=Outputdefinition8Enum();
		data[numpy.nonzero(self.cost_functions==1010)]=Outputdefinition10Enum();
		WriteData(fid,'data',data,'enum',InversionCostFunctionsEnum(),'format','DoubleMat','mattype',3);
		WriteData(fid,'data',num_cost_functions,'enum',InversionNumCostFunctionsEnum(),'format','Integer');
		
