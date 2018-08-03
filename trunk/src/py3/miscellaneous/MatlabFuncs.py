def oshostname():
	import socket

	return socket.gethostname()

def ispc():
	import platform

	if 'Windows' in platform.system():
		return True
	else:
		return False

def ismac():
	import platform

	if 'Darwin' in platform.system():
		return True
	else:
		return False

def strcmp(s1,s2):

	if s1 == s2:
		return True
	else:
		return False

def strncmp(s1,s2,n):

	if s1[0:n] == s2[0:n]:
		return True
	else:
		return False

def strcmpi(s1,s2):

	if s1.lower() == s2.lower():
		return True
	else:
		return False

def strncmpi(s1,s2,n):

	if s1.lower()[0:n] == s2.lower()[0:n]:
		return True
	else:
		return False

def ismember(a,s):
	import numpy

	if not isinstance(s,(tuple,list,dict,numpy.ndarray)):
		s=[s]

	if not isinstance(a,(tuple,list,dict,numpy.ndarray)):
		a=[a]

	if not isinstance(a,numpy.ndarray):
		b=[item in s for item in a]

	else:
		if not isinstance(s,numpy.ndarray):
			b=numpy.empty_like(a)
			for i,item in enumerate(a.flat):
				b.flat[i]=item in s
		else:
			b=numpy.in1d(a.flat,s.flat).reshape(a.shape)

	return b

def det(a):
	import numpy

	if   a.shape==(1,):
		return a[0]
	elif a.shape==(1,1):
		return a[0,0]
	elif a.shape==(2,2):
		return a[0,0]*a[1,1]-a[0,1]*a[1,0]
	else:
		raise TypeError("MatlabFunc.det only implemented for shape (2, 2), not for shape %s." % str(a.shape))

def sparse(ivec,jvec,svec,m=0,n=0,nzmax=0):
	import numpy

	if not m:
		m=numpy.max(ivec)
	if not n:
		n=numpy.max(jvec)

	a=numpy.zeros((m,n))

	for i,j,s in zip(ivec.reshape(-1,order='F'),jvec.reshape(-1,order='F'),svec.reshape(-1,order='F')):
		a[i-1,j-1]+=s

	return a

def heaviside(x):
	import numpy

	y=numpy.zeros_like(x)
	y[numpy.nonzero(x> 0.)]=1.
	y[numpy.nonzero(x==0.)]=0.5

	return y

