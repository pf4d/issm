#! /usr/bin/env python

import os
import sys
import numpy
import math
import struct
import argparse

def BinRead(filin,filout='',verbose=0): #{{{

	print "reading binary file."
	f=open(filin,'rb')

	if filout:
		sys.stdout=open(filout,'w')

	while True:
		try:
			#Step 1: read size of record name
			recordnamesize=struct.unpack('i',f.read(struct.calcsize('i')))[0]
		except struct.error as e:
			print "probable EOF: %s" % e
			break

		print "============================================================================"
		if verbose>2:
			print "\nrecordnamesize = \"%d\"" % (recordnamesize)
		recordname=struct.unpack('%ds' % recordnamesize,f.read(recordnamesize))[0]
		print "field: %s" % recordname

		#Step 2: read the data itself.
		#first read length of record
		reclen=struct.unpack('i',f.read(struct.calcsize('i')))[0]
		if verbose>1:
			print "reclen = %d" % reclen

		#read data code: 
		code=struct.unpack('i',f.read(struct.calcsize('i')))[0]
		#print "code = %d (%s)" % (code,CodeToFormat(code))
		print "Format = %s" % CodeToFormat(code)

		if   code == FormatToCode('Boolean'):
#			bval=struct.unpack('b',f.read(reclen-struct.calcsize('i')))[0]
			bval=struct.unpack('i',f.read(reclen-struct.calcsize('i')))[0]
			print "value = %d" % bval

		elif code == FormatToCode('Integer'):
			ival=struct.unpack('i',f.read(reclen-struct.calcsize('i')))[0]
			print "value = %d" % ival

		elif code == FormatToCode('Double'):
			dval=struct.unpack('d',f.read(reclen-struct.calcsize('i')))[0]
			print "value = %f" % dval

		elif code == FormatToCode('String'):
			strlen=struct.unpack('i',f.read(struct.calcsize('i')))[0]
			if verbose>1:
				print "strlen = %d" % strlen
			sval=struct.unpack('%ds' % strlen,f.read(strlen))[0]
			print "value = '%s'" % sval

		elif code == FormatToCode('BooleanMat'):
			#read matrix type: 
			mattype=struct.unpack('i',f.read(struct.calcsize('i')))[0]
			print "mattype = %d" % mattype

			#now read matrix
			s=[0,0]
			s[0]=struct.unpack('i',f.read(struct.calcsize('i')))[0]
			s[1]=struct.unpack('i',f.read(struct.calcsize('i')))[0]
			print "size = [%dx%d]" % (s[0],s[1])
			data=numpy.zeros((s[0],s[1]))
			for i in xrange(s[0]):
				for j in xrange(s[1]):
					data[i][j]=struct.unpack('d',f.read(struct.calcsize('d')))[0]    #get to the "c" convention, hence the transpose
					if verbose>2: print "data[%d,%d] = %f" % (i,j,data[i][j])

		elif code == FormatToCode('IntMat'):
			#read matrix type: 
			mattype=struct.unpack('i',f.read(struct.calcsize('i')))[0]
			print "mattype = %d" % mattype

			#now read matrix
			s=[0,0]
			s[0]=struct.unpack('i',f.read(struct.calcsize('i')))[0]
			s[1]=struct.unpack('i',f.read(struct.calcsize('i')))[0]
			print "size = [%dx%d]" % (s[0],s[1])
			data=numpy.zeros((s[0],s[1]))
			for i in xrange(s[0]):
				for j in xrange(s[1]):
					data[i][j]=struct.unpack('d',f.read(struct.calcsize('d')))[0]    #get to the "c" convention, hence the transpose
					if verbose>2: print "data[%d,%d] = %f" % (i,j,data[i][j])

		elif code == FormatToCode('DoubleMat'):
			#read matrix type: 
			mattype=struct.unpack('i',f.read(struct.calcsize('i')))[0]
			print "mattype = %d" % mattype

			#now read matrix
			s=[0,0]
			s[0]=struct.unpack('i',f.read(struct.calcsize('i')))[0]
			s[1]=struct.unpack('i',f.read(struct.calcsize('i')))[0]
			print "size = [%dx%d]" % (s[0],s[1])
			data=numpy.zeros((s[0],s[1]))
			for i in xrange(s[0]):
				for j in xrange(s[1]):
					data[i][j]=struct.unpack('d',f.read(struct.calcsize('d')))[0]    #get to the "c" convention, hence the transpose
					if verbose>2: print "data[%d,%d] = %f" % (i,j,data[i][j])

		elif code == FormatToCode('MatArray'):
			f.seek(reclen-4,1)
			print "skipping %d bytes for code %d." % (code, reclen-4)

		elif code == FormatToCode('StringArray'):
			f.seek(reclen-4,1)
			print "skipping %d bytes for code %d." % (code, reclen-4)

		else:
			raise TypeError('BinRead error message: data type: %d not supported yet! (%s)' % (code,recordname))

	f.close()
#}}}
def FormatToCode(format): # {{{
	"""
	This routine takes the format string, and hardcodes it into an integer, which 
	is passed along the record, in order to identify the nature of the dataset being 
	sent.
	"""

	if format=='Boolean':
		code=1
	elif format=='Integer':
		code=2
	elif format=='Double':
		code=3
	elif format=='String':
		code=4
	elif format=='BooleanMat':
		code=5
	elif format=='IntMat':
		code=6
	elif format=='DoubleMat':
		code=7
	elif format=='MatArray':
		code=8
	elif format=='StringArray':
		code=9
	else:
		raise InputError('FormatToCode error message: data type not supported yet!')

	return code
# }}}
def CodeToFormat(code): # {{{
	"""
	This routine takes the format string, and hardcodes it into an integer, which 
	is passed along the record, in order to identify the nature of the dataset being 
	sent.
	"""

	if code==1:
		format='Boolean'
	elif code==2:
		format='Integer'
	elif code==3:
		format='Double'
	elif code==4:
		format='String'
	elif code==5:
		format='BooleanMat'
	elif code==6:
		format='IntMat'
	elif code==7:
		format='DoubleMat'
	elif code==8:
		format='MatArray'
	elif code==9:
		format='StringArray'
	else:
		raise TypeError('FormatToCode error message: code %d not supported yet!' %code)

	return format
# }}}

if __name__ == '__main__': #{{{
	if 'PYTHONSTARTUP' in os.environ:
		PYTHONSTARTUP=os.environ['PYTHONSTARTUP']
		print 'PYTHONSTARTUP =',PYTHONSTARTUP
		if os.path.exists(PYTHONSTARTUP):
			try:
				execfile(PYTHONSTARTUP)
			except Exception as e:
				print "PYTHONSTARTUP error: ",e
		else:
			print "PYTHONSTARTUP file '%s' does not exist." % PYTHONSTARTUP

	parser = argparse.ArgumentParser(description='BinRead - function to read binary input file.')
	parser.add_argument('-f','--filin', help='name of binary input file', default='')
	parser.add_argument('-o','--filout', help='optional name of text output file', default='')
	parser.add_argument('-v','--verbose', help='optional level of output', default=0)
	args = parser.parse_args()

	BinRead(args.filin, args.filout,args.verbose)
#}}}
