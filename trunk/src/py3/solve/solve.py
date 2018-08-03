import datetime
import os
import shutil
from pairoptions import pairoptions
from EnumDefinitions import *
from EnumToString import EnumToString
from ismodelselfconsistent import ismodelselfconsistent
from marshall import marshall
from waitonlock import waitonlock
from loadresultsfromcluster import loadresultsfromcluster
import MatlabFuncs as m

def solve(md,solutionenum,**kwargs):
	"""
	SOLVE - apply solution sequence for this model
 
	   Usage:
	      md=solve(md,solutionenum,varargin)
	      where varargin is a list of paired arguments of string OR enums
 
	   solution types available comprise:
	      - StressbalanceSolutionEnum
	      - MasstransportSolutionEnum
	      - ThermalSolutionEnum
	      - SteadystateSolutionEnum
	      - TransientSolutionEnum
	      - BalancethicknessSolutionEnum
	      - BedSlopeSolutionEnum
	      - SurfaceSlopeSolutionEnum
	      - HydrologySolutionEnum
	      - FlaimSolutionEnum
 
	   extra options:
	      - loadonly : does not solve. only load results
		  - checkconsistency : 'yes' or 'no' (default is 'yes'), ensures checks on consistency of model
		  - restart: 'directory name (relative to the execution directory) where the restart file is located.
 
	   Examples:
	      md=solve(md,StressbalanceSolutionEnum);
	"""

	#recover and process solve options
	if EnumToString(solutionenum)[0][-8:] != 'Solution':
		raise ValueError("solutionenum '%s' not supported!" % EnumToString(solutionenum)[0])
	options=pairoptions(solutionenum=solutionenum,**kwargs)

	#recover some fields
	md.private.solution=solutionenum
	cluster=md.cluster

	#check model consistency
	if m.strcmpi(options.getfieldvalue('checkconsistency','yes'),'yes'):
		print("checking model consistency")
		if solutionenum == FlaimSolutionEnum():
			md.private.isconsistent=True
			md.mesh.checkconsistency(md,solutionenum)
			md.flaim.checkconsistency(md,solutionenum)
			if not md.private.isconsistent:
				raise RuntimeError("Model not consistent, see messages above.")
		else:
			ismodelselfconsistent(md)

	#First, build a runtime name that is unique
	restart=options.getfieldvalue('restart','')
	if restart == 1:
		pass #do nothing
	else:
		if restart:
			md.private.runtimename=restart
		else:
			if options.getfieldvalue('runtimename',True):
				c=datetime.datetime.now()
				md.private.runtimename="%s-%02i-%02i-%04i-%02i-%02i-%02i-%i" % (md.miscellaneous.name,c.month,c.day,c.year,c.hour,c.minute,c.second,os.getpid())
			else:
				md.private.runtimename=md.miscellaneous.name 

	#if running qmu analysis, some preprocessing of dakota files using models
	#fields needs to be carried out. 
	if md.qmu.isdakota:
		md=preqmu(md,options)

	#flaim analysis
	if solutionenum == FlaimSolutionEnum():
		md=flaim_sol(md,options)
		[md.private.solution]=EnumToString(solutionenum)
		return md

	#Do we load results only?
	if options.getfieldvalue('loadonly',False):
		md=loadresultsfromcluster(md)
		return md


	#Write all input files
	marshall(md)                                           # bin file
	md.toolkits.ToolkitsFile(md.miscellaneous.name+'.toolkits')    # toolkits file
	cluster.BuildQueueScript(md.private.runtimename,md.miscellaneous.name,md.private.solution,md.settings.io_gather,md.debug.valgrind,md.debug.gprof,md.qmu.isdakota)    # queue file

	#Stop here if batch mode
	if m.strcmpi(options.getfieldvalue('batch','no'),'yes'):
		print('batch mode requested: not launching job interactively')
		print('launch solution sequence on remote cluster by hand')
		return md

	#Upload all required files: 
	modelname = md.miscellaneous.name
	filelist  = [modelname+'.bin ',modelname+'.toolkits ',modelname+'.queue ']
	if md.qmu.isdakota:
		filelist.append(modelname+'.qmu.in')

	if not restart:
		cluster.UploadQueueJob(md.miscellaneous.name,md.private.runtimename,filelist)
	
	#Launch job
	cluster.LaunchQueueJob(md.miscellaneous.name,md.private.runtimename,filelist,restart)

	#wait on lock
	if md.settings.waitonlock>0:
		#we wait for the done file
		islock=waitonlock(md)
		if islock==0:    #no results to be loaded
			print('The results must be loaded manually with md=loadresultsfromcluster(md).')
		else:            #load results
			print('loading results from cluster')
			md=loadresultsfromcluster(md)

	#post processes qmu results if necessary
	if md.qmu.isdakota:
		if not strncmpi(options['keep'],'y',1):
			shutil.rmtree('qmu'+str(os.getpid()))

	return md
