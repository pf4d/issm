#! /usr/bin/env python

def IdToName(id):
	"""
	IDTONAME- return name of test
 
	    Usage:
	       name=IdToName(id)
	"""

	infile  = open('test' + str(id) + '.py','r')
	file_text  = infile.readlines()

	string='#Test Name:'
	name=file_text[0]
	name=name[len(string)+1:-1]
	return name
