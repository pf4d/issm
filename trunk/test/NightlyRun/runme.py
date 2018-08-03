#!/usr/bin/env python
import os
import numpy as np
from traceback import format_exc
from sys import float_info
from glob import glob
from socket import gethostname

def runme(id=None,exclude=None,benchmark='nightly',procedure='check',output='none',rank=1,numprocs=1):
	"""
	RUNME - test deck for ISSM nightly runs
 
	    In a test deck directory (tests/Vertification/NightlyRun for example)
	    The following command will launch all the existing tests:
	    >> runme()
	    To run the tests 101 and 102:
	    >> runme(id=[101,102])
	    etc...
 
	    Available options:
	       'id'            followed by the list of ids requested
	       'exclude'       ids to be excluded from the test
	       'benchmark'     'all' (all of the tests)
                          'nightly' (nightly run/ daily run)
                          'ismip'  : validation of ismip-hom tests
                          'eismint': validation of eismint tests
                          'thermal': validation of thermal tests
                          'mesh'   : validation of mesh tests
                          'adolc'  : validation of adolc tests
                          'slr'   : validation of slr tests

	       'procedure'     'check' : run the test (default)
	                       'update': update the archive
 
	    Usage:
	       runme(varargin)
 
	    Examples:
	       runme()
	       runme(exclude=101)
	       runme(id=102,procedure='update')
	"""

	from parallelrange import parallelrange
	from IdToName import IdToName
	from arch import archread
	from arch import archwrite
	from arch import archdisp

	#Get ISSM_DIR variable
	ISSM_DIR=os.environ['ISSM_DIR']

	#Process options
	#GET benchmark {{{
	if not benchmark in ['all','nightly','ismip','eismint','thermal','mesh','validation','tranforcing','adolc','slr','referential']:
		print("runme warning: benchmark '{}' not supported, defaulting to test 'nightly'.".format(benchmark))
		benchmark='nightly'
	# }}}
	#GET procedure {{{
	if not procedure in ['check','update']:
		print("runme warning: procedure '{}' not supported, defaulting to test 'check'.".format(procedure))
		procedure='check'
	# }}}
	#GET output {{{
	if not output in ['nightly','none']:
		print("runme warning: output '{}' not supported, defaulting to test 'none'.".format(output))
		output='none'
	# }}}
	#GET RANK and NUMPROCS for multithreaded runs {{{
	if (numprocs<rank):
		numprocs=1
	# }}}
	#GET ids  {{{
	flist=glob('test*.py')    #File name must start with 'test' and must end by '.py' and must be different than 'test.py'
	list_ids=[int(file[4:-3]) for file in flist if not file == 'test.py']    #Keep test id only (skip 'test' and '.py')
	#print 'list_ids =',list_ids

	i1,i2=parallelrange(rank,numprocs,len(list_ids))    #Get tests for this cpu only
	list_ids=list_ids[i1:i2+1]
	#print 'list_ids after parallelrange =',list_ids
	
	if id:
		if isinstance(id,list):
			test_ids=id
		else:
			test_ids=[id]
		test_ids=set(test_ids).intersection(set(list_ids))
	else:
		test_ids=set(list_ids)
		
		#print 'test_ids after list =',test_ids
	# }}}
	#GET exclude {{{
	if exclude:
		if isinstance(exclude,list):
			exclude_ids=exclude
		else:
			exclude_ids=[exclude]
		test_ids=test_ids.difference(set(exclude_ids))
#	print 'test_ids after exclude =',test_ids
	# }}}
	#Process Ids according to benchmarks {{{
	if benchmark=='nightly':
		test_ids=test_ids.intersection(set(range(1,1000)))
	elif benchmark=='validation':
		test_ids=test_ids.intersection(set(range(1001,2000)))
	elif benchmark=='ismip':
		test_ids=test_ids.intersection(set(range(1101,1200)))
	elif benchmark=='eismint':
		test_ids=test_ids.intersection(set(range(1201,1300)))
	elif benchmark=='thermal':
		test_ids=test_ids.intersection(set(range(1301,1400)))
	elif benchmark=='mesh':
		test_ids=test_ids.intersection(set(range(1401,1500)))
	elif benchmark=='tranforcing':
		test_ids=test_ids.intersection(set(range(1501,1503)))
	elif benchmark=='referential':
		test_ids=test_ids.intersection(set(range(1601,1603)))
	elif benchmark=='slr':
		test_ids=test_ids.intersection(set(range(2001,2500)))
	elif benchmark=='adolc':
		test_ids=test_ids.intersection(set(range(3001,3200)))
	#print 'test_ids after benchmark =',test_ids
	test_ids=list(test_ids)
	test_ids.sort()
	#print 'test_ids after sort =',test_ids
	# }}}

	#Loop over tests and launch sequence
	root=os.getcwd()
	for id in test_ids:
		print "----------------starting:%i-----------------------" % id
		try:

			#Execute test
			os.chdir(root)
			id_string=IdToName(id)
			execfile('test'+str(id)+'.py',globals())

			#UPDATE ARCHIVE?
			archive_name='Archive'+str(id)
			if procedure=='update':
				archive_file=os.path.join('..','Archives',archive_name+'.arch')
				if os.path.isfile(archive_file):
					os.remove(archive_file)
				for k,fieldname in enumerate(field_names):
					field=np.array(field_values[k],dtype=float)
					if len(field.shape) == 1:
						if np.size(field):
							field=field.reshape(np.size(field),1)
						else:
							field=field.reshape(0,0)
					elif len(field.shape) == 0:
						field=field.reshape(1,1)
					# Matlab uses base 1, so use base 1 in labels
					archwrite(archive_file,archive_name+'_field'+str(k+1),field)
				print "File '%s' saved.\n" % os.path.join('..','Archives',archive_name+'.arch')

			#ELSE: CHECK TEST
			else:

				#load archive
				if os.path.exists(os.path.join('..','Archives',archive_name+'.arch')):
					archive_file=os.path.join('..','Archives',archive_name+'.arch')
				else:
					raise IOError("Archive file '"+os.path.join('..','Archives',archive_name+'.arch')+"' does not exist.")

				for k,fieldname in enumerate(field_names):
					try:
						#Get field and tolerance
						field=np.array(field_values[k])
						if len(field.shape) == 1:
							if np.size(field):
								field=field.reshape(np.size(field),1)
							else:
								field=field.reshape(0,0)
						tolerance=field_tolerances[k]

						#compare to archive
						# Matlab uses base 1, so use base 1 in labels
						archive=np.array(archread(archive_file,archive_name+'_field'+str(k+1)))
						if archive == None:
							raise NameError("Field name '"+archive_name+'_field'+str(k+1)+"' does not exist in archive file.")
						error_diff=np.amax(np.abs(archive-field),axis=0)/(np.amax(np.abs(archive),axis=0)+float_info.epsilon)

						#disp test result
						if (np.any(error_diff>tolerance) or np.isnan(error_diff)):
							print('ERROR   difference: {} > {} test id: {} test name: {} field: {}'.format(error_diff,tolerance,id,id_string,fieldname))
						else:
							print('SUCCESS difference: {} < {} test id: {} test name: {} field: {}'.format(error_diff,tolerance,id,id_string,fieldname))

					except Exception as message:

						#something went wrong, print failure message:
						print format_exc()
						directory=os.getcwd().split('/')    #  not used?
						if output=='nightly':
							fid=open(os.path.join(ISSM_DIR,'nightlylog','pythonerror.log'), 'a')
							fid.write('%s' % message)
							fid.write('\n------------------------------------------------------------------\n')
							fid.close()
							print('FAILURE difference: N/A test id: {} test name: {} field: {}'.format(id,id_string,fieldname))
						else:
							print('FAILURE difference: N/A test id: {} test name: {} field: {}'.format(id,id_string,fieldname))
							raise RuntimeError(message)


		except Exception as message:

			#something went wrong, print failure message:
			print format_exc()
			directory=os.getcwd().split('/')    #  not used?
			if output=='nightly':
				fid=open(os.path.join(ISSM_DIR,'nightlylog','pythonerror.log'), 'a')
				fid.write('%s' % message)
				fid.write('\n------------------------------------------------------------------\n')
				fid.close()
				print('FAILURE difference: N/A test id: {} test name: {} field: {}'.format(id,id_string,'N/A'))
			else:
				print('FAILURE difference: N/A test id: {} test name: {} field: {}'.format(id,id_string,'N/A'))
				raise RuntimeError(message)

		print "----------------finished:%i-----------------------" % id
	return

import argparse
if __name__ == '__main__':
	if 'PYTHONSTARTUP' in os.environ:
		PYTHONSTARTUP=os.environ['PYTHONSTARTUP']
		#print 'PYTHONSTARTUP =',PYTHONSTARTUP
		if os.path.exists(PYTHONSTARTUP):
			try:
				execfile(PYTHONSTARTUP)
			except Exception as e:
				print "PYTHONSTARTUP error: ",e
		else:
			print("PYTHONSTARTUP file '{}' does not exist.".format(PYTHONSTARTUP))

	parser = argparse.ArgumentParser(description='RUNME - test deck for ISSM nightly runs')
	parser.add_argument('-i','--id', nargs='*', type=int, help='followed by the list of ids requested', default=[])
	parser.add_argument('-e','--exclude', nargs='+', type=int, help='ids to be excluded from the test', default=[])
	parser.add_argument('-b','--benchmark', help='nightly/ismip/eismint/thermal/mesh/...', default='nightly')
	parser.add_argument('-p','--procedure', help='check/update', default='check')
	parser.add_argument('-o','--output', help='nightly/daily/none', default='none')
	parser.add_argument('-r','--rank', type=int, help='rank', default=1)
	parser.add_argument('-n','--numprocs', type=int, help='numprocs', default=1)
	args = parser.parse_args()

	md = runme(args.id, args.exclude, args.benchmark, args.procedure, args.output, args.rank, args.numprocs)

	if args.output=='nightly':
		print "PYTHONEXITEDCORRECTLY"

	exit(md)
