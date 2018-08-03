try:
	import pylab as p
	import matplotlib.pyplot as plt
except ImportError:
	print("could not import pylab, matplotlib has not been installed, no plotting capabilities enabled")

from checkplotoptions import checkplotoptions
from plot_mesh import plot_mesh
from processmesh import processmesh
from processdata import processdata
from plot_unit import plot_unit
from applyoptions import applyoptions

try:
	from osgeo import gdal
	overlaysupport=True
except ImportError:
	print('osgeo/gdal for python not installed, overlay plots are not enabled')
	overlaysupport=False

if overlaysupport:
	from plot_overlay import plot_overlay

def plot_manager(md,options,fig,ax):
	'''
	PLOT_MANAGER - distribute the plots called by plotmodel

	'fig' is a handle to the figure instance created by plotmodel.

	'ax' is a handle to the axes instance created by plotmodel.  This is
	currently generated using matplotlib's AxesGrid toolkit.

	Usage:
		plot_manager(md,options,fig,ax);

	See also: PLOTMODEL, PLOT_UNIT
	'''

	#parse options and get a structure of options
	options=checkplotoptions(md,options)

	#get data to be plotted
	data=options.getfieldvalue('data');

	#initialize plot handle variable
	#handle=None

	# initialize subplot
	#p.subplot(nrows,ncols,i,aspect='equal')

	##basemap plot
	#if options.exist('basemap'):
	#	plot_basemap(md,data,options,nrows,ncols,i)

	#overlay plot
	if options.exist('overlay') and overlaysupport:
		plot_overlay(md,data,options,ax)
		options.addfielddefault('alpha',0.5)
		options.addfielddefault('xlim',[min(md.mesh.x),max(md.mesh.x)])
		options.addfielddefault('ylim',[min(md.mesh.y),max(md.mesh.y)])

	#figure out if this is a special plot
	if isinstance(data,str):

		# convert string to lower case for a case-insensitive comparison
		if data.lower()=='mesh': 
			plot_mesh(md,options,ax)
			applyoptions(md,[],options,fig,ax)
			fig.delaxes(fig.axes[1]) # hack to remove colorbar after the fact
			return
		elif data.lower()=='none':
			print('no data provided to plot (TODO: write plot_none.py)')
			applyoptions(md,[],options,fig,ax)
			return
		else:
			print(("WARNING: '%s' is not implemented or is not a valid string for option 'data'" % data))

	#elif data in vars(md):
	#else:
		#print "'data' not a string, plotting model properties yet to be implemented..."

	#Gridded plot

	#Section plot

	#Profile plot

	#process data and model
	x,y,z,elements,is2d,isplanet=processmesh(md,data,options)
	data2,datatype=processdata(md,data,options)

	#standard plot
	#if not handle:
	#	p.subplot(nrows,ncols,i,aspect='equal')

	#plot unit
	plot_unit(x,y,z,elements,data2,is2d,isplanet,datatype,options,ax)

	#apply all options
	applyoptions(md,data2,options,fig,ax)
	
	#ground overlay on kml plot_unit
