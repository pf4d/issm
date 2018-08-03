#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

import os,sys
import urllib

#Check inputs
if(len(sys.argv)!=3): raise NameError('usage: ./DownloadExternalPackage.py URL localfile')

url=sys.argv[1];
localFile=sys.argv[2]

#Remove file if it already exists
if os.path.exists(localFile):
	print "File "+ localFile +" already exists and will not be downloaded..."
	sys.exit()

#Try to download from url
httpfail=-1
try:
	print "Fetching %s" % localFile
	urllib.urlretrieve(url,localFile)
	httpfail=0
except Exception, e:
	httpfail=1

#Error message in case it failed
if (httpfail):
	failureMessage = '''
===========================================================================
Unable to download package %s from: %s
* If URL specified manually - perhaps there is a typo?
* If your network is disconnected - please reconnect 
* Alternatively, you can download the above URL manually
===========================================================================
''' % (localFile,url)
	raise RuntimeError(failureMessage)
