import numpy as np
from plotoptions import plotoptions

try:
	import pylab as p
	import matplotlib.pyplot as plt
	from mpl_toolkits.axes_grid1 import ImageGrid, AxesGrid
except ImportError:
	print("could not import pylab, matplotlib has not been installed, no plotting capabilities enabled")

from plot_manager import plot_manager
from math import ceil, sqrt

def plotmodel(md,*args):
	'''
	at command prompt, type 'plotdoc' for additional documentation
	'''

	#First process options 
	options=plotoptions(*args)

	#get number of subplots
	subplotwidth=ceil(sqrt(options.numberofplots))
	
	#Get figure number and number of plots
	figurenumber=options.figurenumber
	numberofplots=options.numberofplots

	#get hold
	hold=options.list[0].getfieldvalue('hold',False)

	#if nrows and ncols specified, then bypass
	if options.list[0].exist('nrows'):
		nrows=options.list[0].getfieldvalue('nrows')
		nr=True
	else:
		nrows=np.ceil(numberofplots/subplotwidth)
		nr=False
	
	if options.list[0].exist('ncols'):
		ncols=options.list[0].getfieldvalue('ncols')
		nc=True
	else:
		ncols=int(subplotwidth)
		nc=False
	ncols=int(ncols)
	nrows=int(nrows)
	
	#check that nrows and ncols were given at the same time!
	if not nr==nc:
		raise Exception('error: nrows and ncols need to be specified together, or not at all')
	
	#Go through plots
	if numberofplots:
		
		#if plt.fignum_exists(figurenumber): 
		#	plt.cla()

		#if figsize specified
		if options.list[0].exist('figsize'):
			figsize=options.list[0].getfieldvalue('figsize')
			fig=plt.figure(figurenumber,figsize=(figsize[0],figsize[1]),tight_layout=True)
		else:
			fig=plt.figure(figurenumber,tight_layout=True)
		fig.clf()

		# options needed to define plot grid
		direction=options.list[0].getfieldvalue('direction','row') # row,column
		axes_pad=options.list[0].getfieldvalue('axes_pad',0.25)
		add_all=options.list[0].getfieldvalue('add_all',True) # True,False
		share_all=options.list[0].getfieldvalue('share_all',True) # True,False
		label_mode=options.list[0].getfieldvalue('label_mode','1') # 1,L,all
		cbar_mode=options.list[0].getfieldvalue('cbar_mode','each') # none,single,each
		cbar_location=options.list[0].getfieldvalue('cbar_location','right') # right,top
		cbar_size=options.list[0].getfieldvalue('cbar_size','5%')
		cbar_pad=options.list[0].getfieldvalue('cbar_pad','2.5%') # None or %
		
		axgrid=ImageGrid(fig, 111,
				nrows_ncols=(nrows,ncols),
				direction=direction,
				axes_pad=axes_pad,
				add_all=add_all,
				share_all=share_all,
				label_mode=label_mode,
				cbar_mode=cbar_mode,
				cbar_location=cbar_location,
				cbar_size=cbar_size,
				cbar_pad=cbar_pad
				)

		if cbar_mode=='none':
			for ax in axgrid.cbar_axes: fig._axstack.remove(ax)

		for i in range(numberofplots):
			plot_manager(options.list[i].getfieldvalue('model',md),options.list[i],fig,axgrid[i])

		fig.show()
	else:
		raise Exception('plotmodel error message: no output data found.')
