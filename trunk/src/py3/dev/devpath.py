#!/usr/bin/env python
import os,sys
import warnings

#Recover ISSM_DIR and USERNAME
ISSM_DIR = os.getenv('ISSM_DIRPY3')
USERNAME = os.getenv('USER')
JPL_SVN  = os.getenv('JPL_SVN')
if(ISSM_DIR==None):
	raise NameError('"ISSM_DIR" environment variable is empty! You should define ISSM_DIR in your .cshrc or .bashrc!')
if(JPL_SVN==None):
	warnings.warn('"JPL_SVN" environment variable is empty! add it to your .cshrc or .bashrc if you want to do distant computing')

#Go through src/m and append any directory that contains a *.py file to PATH 
for root,dirs,files in os.walk(ISSM_DIR+ '/src/py3'):
	if '.svn' in dirs:
		dirs.remove('.svn')
	for file in files:
		if file.find(".py") != -1:
			if file.find(".pyc") == -1:
				if root not in sys.path:
					sys.path.append(root)
				
sys.path.append(ISSM_DIR + '/lib')
sys.path.append(ISSM_DIR + '/src/wrappers/python/.libs')
# If using clusters, we need to have the path to the cluster settings directory
if(JPL_SVN!=None):
	if os.path.exists(JPL_SVN + '/usr/' + USERNAME):
		sys.path.append(JPL_SVN + '/usr/' + USERNAME)
	else:
		raise NameError ('cluster settings should be in, '+ JPL_SVN +'/usr/' + USERNAME)

#Manual imports for commonly used functions
#from plotmodel import plotmodel

#c = get_ipython().config
#c.InteractiveShellApp.exec_lines = []
#c.InteractiveShellApp.exec_lines.append('%load_ext autoreload')
#c.InteractiveShellApp.exec_lines.append('%autoreload 2')
#c.InteractiveShellApp.exec_lines.append('print "Warning: disable autoreload in startup.py to improve performance." ')

print("\n  ISSM development path correctly loaded\n\n")
