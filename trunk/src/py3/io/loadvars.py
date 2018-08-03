import shelve
import os.path
from dbm import whichdb

def loadvars(*args):
	"""
	LOADVARS - function to load variables to a file.

	This function loads one or more variables from a file.  The names of the variables
	must be supplied.  If more than one variable is specified, it may be done with
	a list of names or a dictionary of name as keys.  The output type will correspond
	to the input type.  All the variables in the file may be loaded by specifying only
	the file name.

	Usage:
	   a=loadvars('shelve.dat','a')
	   [a,b]=loadvars('shelve.dat',['a','b'])
	   nvdict=loadvars('shelve.dat',{'a':None,'b':None})
	   nvdict=loadvars('shelve.dat')

	"""

	filename=''
	nvdict={}

	if len(args) >= 1 and isinstance(args[0],str):
		filename=args[0]
		if not filename:
			filename='/tmp/shelve.dat'

	else:
		raise TypeError("Missing file name.")

	if   len(args) >= 2 and isinstance(args[1],str):    # (filename,name)
		for name in args[1:]:
			nvdict[name]=None

	elif len(args) == 2 and isinstance(args[1],list):    # (filename,[names])
		for name in args[1]:
			nvdict[name]=None

	elif len(args) == 2 and isinstance(args[1],dict):    # (filename,{names:values})
		nvdict=args[1]

	elif len(args) == 1:    #  (filename)
		pass

	else:
		raise TypeError("Unrecognized input arguments.")

	if whichdb(filename):
		print("Loading variables from file '%s'." % filename)
	else:
		raise IOError("File '%s' not found." % filename)

	my_shelf = shelve.open(filename,'r') # 'r' for read-only

	if nvdict:
		for name in nvdict.keys():
			try:
				nvdict[name] = my_shelf[name]
				print("Variable '%s' loaded." % name)
			except KeyError:
				value = None
				print("Variable '%s' not found." % name)

	else:
		for name in my_shelf.keys():
			nvdict[name] = my_shelf[name]
			print("Variable '%s' loaded." % name)

	my_shelf.close()

	if   len(args) >= 2 and isinstance(args[1],str):    # (value)
		value=[nvdict[name] for name in args[1:]]
		return value

	elif len(args) == 2 and isinstance(args[1],list):    # ([values])
		value=[nvdict[name] for name in args[1]]
		return value

	elif (len(args) == 2 and isinstance(args[1],dict)) or (len(args) == 1):    # ({names:values})
		return nvdict

