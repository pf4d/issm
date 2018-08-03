from EnumDefinitions import *
from EnumToString import EnumToString

def AnalysisConfiguration(solutiontype): #{{{
	"""
	ANALYSISCONFIGURATION - return type of analyses, number of analyses 

		Usage:
			[analyses]=AnalysisConfiguration(solutiontype);
	"""

	if   solutiontype == StressbalanceSolutionEnum():
		analyses=[StressbalanceAnalysisEnum(),StressbalanceVerticalAnalysisEnum(),StressbalanceSIAAnalysisEnum(),L2ProjectionBaseAnalysisEnum()]

	elif solutiontype == SteadystateSolutionEnum():
		analyses=[StressbalanceAnalysisEnum(),StressbalanceVerticalAnalysisEnum(),StressbalanceSIAAnalysisEnum(),L2ProjectionBaseAnalysisEnum(),ThermalAnalysisEnum(),MeltingAnalysisEnum()]

	elif solutiontype == ThermalSolutionEnum():
		analyses=[EnthalpyAnalysisEnum(),ThermalAnalysisEnum(),MeltingAnalysisEnum()]

	elif solutiontype == MasstransportSolutionEnum():
		analyses=[MasstransportAnalysisEnum()]

	elif solutiontype == BalancethicknessSolutionEnum():
		analyses=[BalancethicknessAnalysisEnum()]

	elif solutiontype == SurfaceSlopeSolutionEnum():
		analyses=[L2ProjectionBaseAnalysisEnum()]

	elif solutiontype == BalancevelocitySolutionEnum():
		analyses=[BalancevelocityAnalysisEnum()]

	elif solutiontype == BedSlopeSolutionEnum():
		analyses=[L2ProjectionBaseAnalysisEnum()]

	elif solutiontype == GiaSolutionEnum():
		analyses=[GiaAnalysisEnum()]

	elif solutiontype == TransientSolutionEnum():
		analyses=[StressbalanceAnalysisEnum(),StressbalanceVerticalAnalysisEnum(),StressbalanceSIAAnalysisEnum(),L2ProjectionBaseAnalysisEnum(),ThermalAnalysisEnum(),MeltingAnalysisEnum(),EnthalpyAnalysisEnum(),MasstransportAnalysisEnum()]

	elif solutiontype == FlaimSolutionEnum():
		analyses=[FlaimAnalysisEnum()]

	elif solutiontype == HydrologySolutionEnum():
		analyses=[L2ProjectionBaseAnalysisEnum(),HydrologyShreveAnalysisEnum(),HydrologyDCInefficientAnalysisEnum(),HydrologyDCEfficientAnalysisEnum()]

	elif DamageEvolutionSolutionEnum():
		analyses=[DamageEvolutionAnalysisEnum()]

	else:
		raise TypeError("solution type: '%s' not supported yet!" % EnumToString(solutiontype)[0])

	return analyses
#}}}

def ismodelselfconsistent(md):
	"""
	ISMODELSELFCONSISTENT - check that model forms a closed form solvable problem.

	   Usage:
	      ismodelselfconsistent(md),
	"""

	#initialize consistency as true
	md.private.isconsistent=True

	#Get solution and associated analyses
	solution=md.private.solution
	analyses=AnalysisConfiguration(solution)

	#Go through a model fields, check that it is a class, and call checkconsistency
	fields=vars(md)
#	for field in fields.iterkeys():
	for field in md.properties():

		#Some properties do not need to be checked
		if field in ['results','debug','radaroverlay']:
			continue

		#Check that current field is an object
		if not hasattr(getattr(md,field),'checkconsistency'):
			md.checkmessage("field '%s' is not an object." % field)

		#Check consistency of the object
		exec("md.%s.checkconsistency(md,solution,analyses)" % field)

	#error message if mode is not consistent
	if not md.private.isconsistent:
		raise RuntimeError('Model not consistent, see messages above.')

