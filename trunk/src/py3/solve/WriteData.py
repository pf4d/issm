import numpy
import math
import struct
import pairoptions
import MatlabFuncs as m
from EnumDefinitions import *
from EnumToString import EnumToString

def WriteData(fid,**kwargs):
	"""
	WRITEDATA - write model field in binary file
 
	   Usage:
	      WriteData(fid,varargin)
	"""

	#process options
	options=pairoptions.pairoptions(**kwargs)

	#Get data properties
	if options.exist('object'):
		#This is an object field, construct enum and data
		obj       = options.getfieldvalue('object')
		fieldname = options.getfieldvalue('fieldname')
		classname = options.getfieldvalue('class',str(type(obj)).rsplit('.')[-1].split("'")[0])
		if options.exist('enum'):
			enum = options.getfieldvalue('enum')
		else:
			enum = BuildEnum(classname+'_'+fieldname)
		data      = getattr(obj,fieldname)
	else:
		#No processing required
		data = options.getfieldvalue('data')
		enum = options.getfieldvalue('enum')
	format  = options.getfieldvalue('format')
	mattype = options.getfieldvalue('mattype',0)    #only required for matrices
	timeserieslength = options.getfieldvalue('timeserieslength',-1)

	#Process sparse matrices
#	if issparse(data),
#		data=full(data);
#	end

	#Scale data if necesarry
	if options.exist('scale'):
		scale = options.getfieldvalue('scale')
		if numpy.size(data) > 1 :
			if numpy.size(data,0)==timeserieslength:
				data=numpy.array(data)
				data[0:-1,:] = scale*data[0:-1,:]
			else:
				data  = scale*data
		else:
			data  = scale*data
	if numpy.size(data) > 1 :
		if numpy.size(data,0)==timeserieslength:
			yts=365.0*24.0*3600.0
			data[-1,:] = yts*data[-1,:]

	#Step 1: write the enum to identify this record uniquely
	fid.write(struct.pack('i',enum)) 

	#Step 2: write the data itself.
	if   m.strcmpi(format,'Boolean'):    # {{{
#		if len(data) !=1:
#			raise ValueError('field %s cannot be marshalled as it has more than one element!' % EnumToString(enum)[0])

		#first write length of record
		fid.write(struct.pack('i',4+4))  #1 bool (disguised as an int)+code

		#write data code: 
		fid.write(struct.pack('i',FormatToCode(format))) 

		#now write integer
		fid.write(struct.pack('i',int(data)))  #send an int, not easy to send a bool
		# }}}

	elif m.strcmpi(format,'Integer'):    # {{{
#		if len(data) !=1:
#			raise ValueError('field %s cannot be marshalled as it has more than one element!' % EnumToString(enum)[0])

		#first write length of record
		fid.write(struct.pack('i',4+4))  #1 integer + code

		#write data code: 
		fid.write(struct.pack('i',FormatToCode(format))) 

		#now write integer
		fid.write(struct.pack('i',data)) 
		# }}}

	elif m.strcmpi(format,'Double'):    # {{{
#		if len(data) !=1:
#			raise ValueError('field %s cannot be marshalled as it has more than one element!' % EnumToString(enum)[0])

		#first write length of record
		fid.write(struct.pack('i',8+4))  #1 double+code

		#write data code: 
		fid.write(struct.pack('i',FormatToCode(format))) 

		#now write double
		fid.write(struct.pack('d',data)) 
		# }}}

	elif m.strcmpi(format,'String'):    # {{{
		#first write length of record
		fid.write(struct.pack('i',len(data)+4+4))  #string + string size + code

		#write data code: 
		fid.write(struct.pack('i',FormatToCode(format))) 

		#now write string
		fid.write(struct.pack('i',len(data))) 
		fid.write(struct.pack('%ds' % len(data),data)) 
		# }}}

	elif m.strcmpi(format,'BooleanMat'):    # {{{

		if   isinstance(data,bool):
			data=numpy.array([data])
		elif isinstance(data,(list,tuple)):
			data=numpy.array(data).reshape(-1,1)
		if numpy.ndim(data) == 1:
			if numpy.size(data):
				data=data.reshape(numpy.size(data),1)
			else:
				data=data.reshape(0,0)

		#Get size
		s=data.shape
		#if matrix = NaN, then do not write anything
		if s[0]==1 and s[1]==1 and math.isnan(data[0][0]):
			s=(0,0)

		#first write length of record
		fid.write(struct.pack('i',4+4+8*s[0]*s[1]+4+4))    #2 integers (32 bits) + the double matrix + code + matrix type

		#write data code and matrix type: 
		fid.write(struct.pack('i',FormatToCode(format))) 
		fid.write(struct.pack('i',mattype))

		#now write matrix
		fid.write(struct.pack('i',s[0])) 
		fid.write(struct.pack('i',s[1])) 
		for i in range(s[0]):
			for j in range(s[1]):
				fid.write(struct.pack('d',float(data[i][j])))    #get to the "c" convention, hence the transpose
		# }}}

	elif m.strcmpi(format,'IntMat'):    # {{{

		if   isinstance(data,int):
			data=numpy.array([data])
		elif isinstance(data,(list,tuple)):
			data=numpy.array(data).reshape(-1,1)
		if numpy.ndim(data) == 1:
			if numpy.size(data):
				data=data.reshape(numpy.size(data),1)
			else:
				data=data.reshape(0,0)

		#Get size
		s=data.shape
		#if matrix = NaN, then do not write anything
		if s[0]==1 and s[1]==1 and math.isnan(data[0][0]):
			s=(0,0)

		#first write length of record
		fid.write(struct.pack('i',4+4+8*s[0]*s[1]+4+4))    #2 integers (32 bits) + the double matrix + code + matrix type

		#write data code and matrix type: 
		fid.write(struct.pack('i',FormatToCode(format))) 
		fid.write(struct.pack('i',mattype))

		#now write matrix
		fid.write(struct.pack('i',s[0])) 
		fid.write(struct.pack('i',s[1])) 
		for i in range(s[0]):
			for j in range(s[1]):
				fid.write(struct.pack('d',float(data[i][j])))    #get to the "c" convention, hence the transpose
		# }}}

	elif m.strcmpi(format,'DoubleMat'):    # {{{

		if   isinstance(data,(bool,int,float)):
			data=numpy.array([data])
		elif isinstance(data,(list,tuple)):
			data=numpy.array(data).reshape(-1,1)
		if numpy.ndim(data) == 1:
			if numpy.size(data):
				data=data.reshape(numpy.size(data),1)
			else:
				data=data.reshape(0,0)

		#Get size
		s=data.shape
		#if matrix = NaN, then do not write anything
		if s[0]==1 and s[1]==1 and math.isnan(data[0][0]):
			s=(0,0)

		#first write length of record
		recordlength=4+4+8*s[0]*s[1]+4+4; #2 integers (32 bits) + the double matrix + code + matrix type
		if recordlength > 2**31 :
			raise ValueError('field %s cannot be marshalled because it is larger than 4^31 bytes!' % EnumToString(enum)[0])

		fid.write(struct.pack('i',recordlength))  #2 integers (32 bits) + the double matrix + code + matrix type

		#write data code and matrix type: 
		fid.write(struct.pack('i',FormatToCode(format))) 
		fid.write(struct.pack('i',mattype))

		#now write matrix
		fid.write(struct.pack('i',s[0])) 
		fid.write(struct.pack('i',s[1])) 
		for i in range(s[0]):
			for j in range(s[1]):
				fid.write(struct.pack('d',float(data[i][j])))    #get to the "c" convention, hence the transpose
		# }}}

	elif m.strcmpi(format,'MatArray'):    # {{{

		#first get length of record
		recordlength=4+4    #number of records + code
		for matrix in data:
			if   isinstance(matrix,(bool,int,float)):
				matrix=numpy.array([matrix])
			elif isinstance(matrix,(list,tuple)):
				matrix=numpy.array(matrix).reshape(-1,1)
			if numpy.ndim(matrix) == 1:
				if numpy.size(matrix):
					matrix=matrix.reshape(numpy.size(matrix),1)
				else:
					matrix=matrix.reshape(0,0)

			s=matrix.shape
			recordlength+=4*2+s[0]*s[1]*8    #row and col of matrix + matrix of doubles

		#write length of record
		fid.write(struct.pack('i',recordlength)) 

		#write data code: 
		fid.write(struct.pack('i',FormatToCode(format))) 

		#write data, first number of records
		fid.write(struct.pack('i',len(data))) 

		#write each matrix: 
		for matrix in data:
			if   isinstance(matrix,(bool,int,float)):
				matrix=numpy.array([matrix])
			elif isinstance(matrix,(list,tuple)):
				matrix=numpy.array(matrix).reshape(-1,1)
			if numpy.ndim(matrix) == 1:
				matrix=matrix.reshape(numpy.size(matrix),1)

			s=matrix.shape
			fid.write(struct.pack('i',s[0])) 
			fid.write(struct.pack('i',s[1])) 
			for i in range(s[0]):
				for j in range(s[1]):
					fid.write(struct.pack('d',float(matrix[i][j])))
		# }}}

	elif m.strcmpi(format,'StringArray'):    # {{{

		#first get length of record
		recordlength=4+4    #for length of array + code
		for string in data:
			recordlength+=4+len(string)    #for each string

		#write length of record
		fid.write(struct.pack('i',recordlength)) 

		#write data code: 
		fid.write(struct.pack('i',FormatToCode(format))) 

		#now write length of string array
		fid.write(struct.pack('i',len(data))) 

		#now write the strings
		for string in data:
			fid.write(struct.pack('i',len(string))) 
			fid.write(struct.pack('%ds' % len(string),string)) 
		# }}}

	else:    # {{{
		raise TypeError('WriteData error message: data type: %d not supported yet! (%s)' % (format,EnumToString(enum)[0]))
	# }}}

def BuildEnum(string): # {{{
	"""
	BUILDENUM - build enum out of string
 
    Usage:
       enum=BuildEnum(string)
	"""

	if '_' in string:
		substrs=string.split('_')
		string=''
		for substr in substrs:
			string+=substr[0].upper()+substr[1:]
	else:
		#take first letter of string and make it uppercase: 
		string=string[0].upper()+string[1:]

	#Get Enum
	enum=StringToEnum(string)[0]

	return enum
# }}}

def FormatToCode(format): # {{{
	"""
	This routine takes the format string, and hardcodes it into an integer, which 
	is passed along the record, in order to identify the nature of the dataset being 
	sent.
	"""

	if   m.strcmpi(format,'Boolean'):
		code=1
	elif m.strcmpi(format,'Integer'):
		code=2
	elif m.strcmpi(format,'Double'):
		code=3
	elif m.strcmpi(format,'String'):
		code=4
	elif m.strcmpi(format,'BooleanMat'):
		code=5
	elif m.strcmpi(format,'IntMat'):
		code=6
	elif m.strcmpi(format,'DoubleMat'):
		code=7
	elif m.strcmpi(format,'MatArray'):
		code=8
	elif m.strcmpi(format,'StringArray'):
		code=9
	else:
		raise InputError('FormatToCode error message: data type not supported yet!')

	return code
# }}}

