import numpy

def ComputeMetric(hessian,scale,epsilon,hmin,hmax,pos):
	"""
	COMPUTEMETRIC - compute metric from an Hessian

	   Usage:
	      metric=ComputeMetric(hessian,scale,epsilon,hmin,hmax,pos)
	      pos is contains the positions where the metric is wished to be maximized (water?)

	   Example:
	      metric=ComputeMetric(hessian,2/9,10^-1,100,10^5,[])
	"""

	#first, find the eigen values of each line of H=[hessian(i,1) hessian(i,2); hessian(i,2) hessian(i,3)]
	a=hessian[:,0]
	b=hessian[:,1]
	d=hessian[:,2]
	lambda1=0.5*((a+d)+numpy.sqrt(4.*b**2+(a-d)**2))
	lambda2=0.5*((a+d)-numpy.sqrt(4.*b**2+(a-d)**2))
	pos1=numpy.nonzero(lambda1==0.)[0]
	pos2=numpy.nonzero(lambda2==0.)[0]
	pos3=numpy.nonzero(numpy.logical_and(b==0.,lambda1==lambda2))[0]

	#Modify the eigen values to control the shape of the elements
	lambda1=numpy.minimum(numpy.maximum(numpy.abs(lambda1)*scale/epsilon,1./hmax**2),1./hmin**2)
	lambda2=numpy.minimum(numpy.maximum(numpy.abs(lambda2)*scale/epsilon,1./hmax**2),1./hmin**2)

	#compute eigen vectors
	norm1=numpy.sqrt(8.*b**2+2.*(d-a)**2+2.*(d-a)*numpy.sqrt((a-d)**2+4.*b**2))
	v1x=2.*b/norm1
	v1y=((d-a)+numpy.sqrt((a-d)**2+4.*b**2))/norm1
	norm2=numpy.sqrt(8.*b**2+2.*(d-a)**2-2.*(d-a)*numpy.sqrt((a-d)**2+4.*b**2))
	v2x=2.*b/norm2
	v2y=((d-a)-numpy.sqrt((a-d)**2+4.*b**2))/norm2

	v1x[pos3]=1.
	v1y[pos3]=0.
	v2x[pos3]=0.
	v2y[pos3]=1.

	#Compute new metric (for each node M=V*Lambda*V^-1)
	metric=numpy.hstack((((v1x*v2y-v1y*v2x)**(-1)*( lambda1*v2y*v1x-lambda2*v1y*v2x)).reshape(-1,1), \
		                 ((v1x*v2y-v1y*v2x)**(-1)*( lambda1*v1y*v2y-lambda2*v1y*v2y)).reshape(-1,1), \
		                 ((v1x*v2y-v1y*v2x)**(-1)*(-lambda1*v2x*v1y+lambda2*v1x*v2y)).reshape(-1,1)))

	#some corrections for 0 eigen values
	metric[pos1,:]=numpy.tile(numpy.array([[1./hmax**2,0.,1./hmax**2]]),(numpy.size(pos1),1))
	metric[pos2,:]=numpy.tile(numpy.array([[1./hmax**2,0.,1./hmax**2]]),(numpy.size(pos2),1))

	#take care of water elements
	metric[pos ,:]=numpy.tile(numpy.array([[1./hmax**2,0.,1./hmax**2]]),(numpy.size(pos ),1))

	#take care of NaNs if any (use Numpy eig in a loop)
	pos=numpy.nonzero(numpy.isnan(metric))[0]
	if numpy.size(pos):
		print((" %i NaN found in the metric. Use Numpy routine..." % numpy.size(pos)))
		for posi in pos:
			H=numpy.array([[hessian[posi,0],hessian[posi,1]],[hessian[posi,1],hessian[posi,2]]])
			[v,u]=numpy.linalg.eig(H)
			v=numpy.diag(v)
			lambda1=v[0,0]
			lambda2=v[1,1]
			v[0,0]=numpy.minimum(numpy.maximum(numpy.abs(lambda1)*scale/epsilon,1./hmax**2),1./hmin**2)
			v[1,1]=numpy.minimum(numpy.maximum(numpy.abs(lambda2)*scale/epsilon,1./hmax**2),1./hmin**2)

			metricTria=numpy.dot(numpy.dot(u,v),numpy.linalg.inv(u))
			metric[posi,:]=numpy.array([metricTria[0,0],metricTria[0,1],metricTria[1,1]])

	if numpy.any(numpy.isnan(metric)):
		raise RunTimeError("ComputeMetric error message: NaN in the metric despite our efforts...")

	return metric

