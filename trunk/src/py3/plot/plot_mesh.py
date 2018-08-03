try:
	import pylab as p
except ImportError:
	print("could not import pylab, matplotlib has not been installed, no plotting capabilities enabled")

from processmesh import processmesh
from applyoptions import applyoptions

def plot_mesh(md,options,ax):
	'''
	PLOT_MESH - plot model mesh

		Usage:
			plot_mesh(md,options,nlines,ncols,i)

		See also: PLOTMODEL
	'''

	x,y,z,elements,is2d,isplanet=processmesh(md,[],options)

	if is2d:
		ax.triplot(x,y,elements)
	else:
		print('WARNING: only 2D mesh plot is currently implemented')
	
	#apply options
	options.addfielddefault('title','Mesh')
	options.addfielddefault('colorbar','off')
	options.addfielddefault('ticklabels','on')
