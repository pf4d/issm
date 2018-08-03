import numpy as  np

try:
	import matplotlib as mpl
except ImportError:
	print 'cannot import matplotlib, no plotting capabilities enabled'

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
	'''
	truncate a colormap within normalized limits [0,1]

	cmap - a matplotlib colormap
	minval - minimum value, normalized, of cmap to be returned.
	maxval - maximum value, normalized, of cmap to be returned.
	n - number of levels to use in constructing the new colormap

	Example:
		newcmap=truncate_colormap(oldcmap,minval=0.2,maxval=0.8,n=128)

	'''

	new_cmap = mpl.colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name,
		a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
	
	return new_cmap
