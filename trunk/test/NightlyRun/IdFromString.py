#! /usr/bin/env python
from IdToName import IdToName

def IdFromString(string):
	"""
	IDFROMSTRING - output ids from a given string
 
	    Usage:
	       ids=IdFromString(string)
 
	    Examples:
	       ids=IdFromString('Parallel')
	       ids=IdFromString('79North')
			 ids=IdFromString('*')
	"""

	#Check input
	if not isinstance(string,str):
		raise TypeError('IdFromString error message: input argument is not a string.')

	#Get the dictionary and scan for matches
	idnames=IdToName(0)
	ids=[item[0] for item in idnames.iteritems() if string in item[1]]

	#Return if no test found
	if not ids:
		print "No test matches '%s'." % string
		return ids

	#Display names
	ids.sort()
	print "%d tests match '%s':" % (len(ids),string)
	for id in ids:
		print "   %d : %s" % (id,idnames[id])

	return ids

