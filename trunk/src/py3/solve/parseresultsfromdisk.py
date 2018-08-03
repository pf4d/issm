import struct
import numpy
from collections import OrderedDict
import results as resultsclass
import MatlabFuncs as m

def parseresultsfromdisk(filename,iosplit):
	"""
	PARSERESULTSFROMDISK - ...

	   Usage:
	      results=parseresultsfromdisk(filename,iosplit)
	"""

	if iosplit:
		results=parseresultsfromdiskiosplit(filename)
	else:
		results=parseresultsfromdiskioserial(filename)

	return results

def parseresultsfromdiskioserial(filename):    # {{{
	"""
	PARSERESULTSFROMDISK - ...
	 
	    Usage:
	       results=parseresultsfromdiskioserial(filename)
	"""

	#Open file
	try:
		fid=open(filename,'rb')
	except IOError as e:
		raise IOError("loadresultsfromdisk error message: could not open '%s' for binary reading." % filename)

	#initialize results: 
	results=[]
	results.append(None)

	#Read fields until the end of the file.
	result=ReadData(fid)

	counter=0
	check_nomoresteps=0
	step=result['step']

	while result:

		if check_nomoresteps:
			#check that the new result does not add a step, which would be an error: 
			if result['step']>=1:
				raise TypeError("parsing results for a steady-state core, which incorporates transient results!")

		#Check step, increase counter if this is a new step
		if(step!=result['step'] and result['step']>1):
			counter = counter + 1
			step    = result['step']

		#Add result
		if result['step']==0:
			#if we have a step = 0, this is a steady state solution, don't expect more steps. 
			index = 0;
			check_nomoresteps=1
	
		elif result['step']==1:
			index = 0
		else:
			index = counter;
	
		if index > len(results)-1:
			for i in range(len(results)-1,index-1):
				results.append(None)
			results.append(resultsclass.results())
		
		elif results[index] is None:
			results[index]=resultsclass.results()

			
		#Get time and step
		if result['step'] != -9999.:
			setattr(results[index],'step',result['step'])
		if result['time'] != -9999.:
			setattr(results[index],'time',result['time']) 
	
		#Add result
		if hasattr(results[index],result['fieldname']) and not m.strcmp(result['fieldname'],'SolutionType'):
			setattr(results[index],result['fieldname'],numpy.vstack((getattr(results[index],result['fieldname']),result['field'])))
		else:
			setattr(results[index],result['fieldname'],result['field'])

		#read next result
		result=ReadData(fid)

	fid.close()

	return results
	# }}}
def parseresultsfromdiskiosplit(filename):    # {{{
	"""
	PARSERESULTSFROMDISKIOSPLIT - ...
	 
	    Usage:
	       results=parseresultsfromdiskiosplit(filename)
	"""

	#Open file
	try:
		fid=open(filename,'rb')
	except IOError as e:
		raise IOError("loadresultsfromdisk error message: could not open '%s' for binary reading." % filename)

	results=[]

	#if we have done split I/O, ie, we have results that are fragmented across patches, 
	#do a first pass, and figure out the structure of results
	result=ReadDataDimensions(fid)
	while result:

		#Get time and step
		if result['step'] > len(results):
			for i in range(len(results),result['step']-1):
				results.append(None)
			results.append(resultsclass.results())
		setattr(results[result['step']-1],'step',result['step'])
		setattr(results[result['step']-1],'time',result['time']) 

		#Add result
		setattr(results[result['step']-1],result['fieldname'],float('NaN'))

		#read next result
		result=ReadDataDimensions(fid)

	#do a second pass, and figure out the size of the patches
	fid.seek(0)    #rewind
	result=ReadDataDimensions(fid)
	while result:

		#read next result
		result=ReadDataDimensions(fid)

	#third pass, this time to read the real information
	fid.seek(0)    #rewind
	result=ReadData(fid)
	while result:

		#Get time and step
		if result['step']> len(results):
			for i in range(len(results),result['step']-1):
				results.append(None)
			results.append(resultsclass.results())
		setattr(results[result['step']-1],'step',result['step'])
		setattr(results[result['step']-1],'time',result['time']) 

		#Add result
		setattr(results[result['step']-1],result['fieldname'],result['field'])

		#read next result
		result=ReadData(fid)

	#close file
	fid.close()

	return results
	# }}}
def ReadData(fid):    # {{{
	"""
	READDATA - ...
	 
	    Usage:
	       field=ReadData(fid)
	"""

	#read field
	try:
		length=struct.unpack('i',fid.read(struct.calcsize('i')))[0]

		fieldname=struct.unpack('%ds' % length,fid.read(length))[0][:-1]
		time=struct.unpack('d',fid.read(struct.calcsize('d')))[0]
		step=struct.unpack('i',fid.read(struct.calcsize('i')))[0]

		type=struct.unpack('i',fid.read(struct.calcsize('i')))[0]
		M=struct.unpack('i',fid.read(struct.calcsize('i')))[0]
		if   type==1:
			field=numpy.array(struct.unpack('%dd' % M,fid.read(M*struct.calcsize('d'))),dtype=float)
		elif type==2:
			field=struct.unpack('%ds' % M,fid.read(M))[0][:-1]
		elif type==3:
			N=struct.unpack('i',fid.read(struct.calcsize('i')))[0]
#			field=transpose(fread(fid,[N M],'double'));
			field=numpy.zeros(shape=(M,N),dtype=float)
			for i in range(M):
				field[i,:]=struct.unpack('%dd' % N,fid.read(N*struct.calcsize('d')))
		else:
			raise TypeError("cannot read data of type %d" % type)

		#Process units here FIXME: this should not be done here!
		yts=365.0*24.0*3600.0
		if m.strcmp(fieldname,'BalancethicknessThickeningRate'):
			field = field*yts
		elif m.strcmp(fieldname,'Time'):
			field = field/yts
		elif m.strcmp(fieldname,'HydrologyWaterVx'):
			field = field*yts
		elif m.strcmp(fieldname,'HydrologyWaterVy'):
			field = field*yts
		elif m.strcmp(fieldname,'Vx'):
			field = field*yts
		elif m.strcmp(fieldname,'Vy'):
			field = field*yts
		elif m.strcmp(fieldname,'Vz'):
			field = field*yts
		elif m.strcmp(fieldname,'Vel'):
			field = field*yts
		elif m.strcmp(fieldname,'BasalforcingsGroundediceMeltingRate'):
			field = field*yts
		elif m.strcmp(fieldname,'TotalSmb'):
			field = field/10.**12.*yts #(GigaTon/year)
		elif m.strcmp(fieldname,'SmbMassBalance'):
			field = field*yts
		elif m.strcmp(fieldname,'CalvingCalvingrate'):
			field = field*yts


		result=OrderedDict()
		result['fieldname']=fieldname
		result['time']=time
		result['step']=step
		result['field']=field

	except struct.error as e:
		result=None

	return result
	# }}}
def ReadDataDimensions(fid):    # {{{
	"""
	READDATADIMENSIONS - read data dimensions, step and time, but not the data itself.
	 
	    Usage:
	       field=ReadDataDimensions(fid)
	"""

	#read field
	try:
		length=struct.unpack('i',fid.read(struct.calcsize('i')))[0]

		fieldname=struct.unpack('%ds' % length,fid.read(length))[0][:-1]
		time=struct.unpack('d',fid.read(struct.calcsize('d')))[0]
		step=struct.unpack('i',fid.read(struct.calcsize('i')))[0]

		type=struct.unpack('i',fid.read(struct.calcsize('i')))[0]
		M=struct.unpack('i',fid.read(struct.calcsize('i')))[0]
		N=1    #default
		if   type==1:
			fid.seek(M*8,1)
		elif type==2:
			fid.seek(M,1)
		elif type==3:
			N=struct.unpack('i',fid.read(struct.calcsize('i')))[0]
			fid.seek(N*M*8,1)
		else:
			raise TypeError("cannot read data of type %d" % type)

		result=OrderedDict()
		result['fieldname']=fieldname
		result['time']=time
		result['step']=step
		result['M']=M
		result['N']=N

	except struct.error as e:
		result=None

	return result
	# }}}
