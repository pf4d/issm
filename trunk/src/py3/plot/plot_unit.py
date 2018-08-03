from cmaptools import truncate_colormap
try:
    import pylab as p
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    print("could not import pylab, matplotlib has not been installed, no plotting capabilities enabled")

def plot_unit(x,y,z,elements,data,is2d,isplanet,datatype,options,ax):
    """
    PLOT_UNIT - unit plot, display data
    
    	Usage:
    		plot_unit(x,y,z,elements,data,is2d,isplanet,datatype,options)
    
    	See also: PLOTMODEL, PLOT_MANAGER
    """
    
    #edgecolor
    edgecolor=options.getfieldvalue('edgecolor','None')
    
    #number of colorlevels for plots
    colorlevels=options.getfieldvalue('colorlevels',128)
    
    alpha=options.getfieldvalue('alpha',1)
    
    #colormap
    # default sequential colormap
    defaultmap=truncate_colormap(mpl.cm.gnuplot2,0.1,0.9,128)
    cmap=options.getfieldvalue('colormap',defaultmap)
    if options.exist('cmap_set_over'):
        over=options.getfieldvalue('cmap_set_over','0.5')
        cmap.set_over(over)
    if options.exist('cmap_set_under'):
        under=options.getfieldvalue('cmap_set_under','0.5')
        cmap.set_under(under)
    
    #normalize colormap if clim/caxis specified
    if options.exist('clim'):
        lims=options.getfieldvalue('clim',[np.amin(data),np.amax(data)])
    elif options.exist('caxis'):
        lims=options.getfieldvalue('caxis',[np.amin(data),np.amax(data)])
    else:
        if np.amin(data)==np.amax(data):
            lims=[np.amin(data)-0.5,np.amax(data)+0.5]
        else:
    	    lims=[np.amin(data),np.amax(data)]
    norm = mpl.colors.Normalize(vmin=lims[0], vmax=lims[1])
    if datatype==1:
       #element plot
        if is2d:
    	    tri=ax.tripcolor(x,y,elements,data,colorlevels,cmap=cmap,norm=norm,alpha=alpha,edgecolors=edgecolor)
    	else:
    	    raise ValueError('plot_unit error: 3D element plot not supported yet')
    	return 
    
    elif datatype==2:
    	#node plot
    	if is2d:
    	    tri=ax.tricontourf(x,y,elements,data,colorlevels,cmap=cmap,norm=norm,alpha=alpha,extend='both')
    	    if edgecolor != 'None':
    	        ax.triplot(x,y,elements,color=edgecolor)
    	else:
    	    raise ValueError('plot_unit error: 3D node plot not supported yet')
    	return
    
    elif datatype==3:
        vx=data[:,0]
        vy=data[:,1]
        #TODO write plot_quiver.py to handle this here
        color=np.sqrt(vx**2+vy**2)
        scale=options.getfieldvalue('scale',1000)
        width=options.getfieldvalue('width',0.005*(np.amax(x)-np.amin(y)))
        headwidth=options.getfieldvalue('headwidth',3)
        headlength=options.getfieldvalue('headlength',5)
        Q=ax.quiver(x,y,vx,vy,color,cmap=cmap,norm=norm,scale=scale,
                width=width,headwidth=headwidth,headlength=headlength)
    	return
    
    elif datatype==4:
    	#P1 patch plot
    	print('plot_unit message: P1 patch plot not implemented yet')
    	return
    
    elif datatype==5:
    	print('plot_unit message: P0 patch plot not implemented yet')
    	return
    
    else:
    	raise ValueError('datatype=%d not supported' % datatype)
    
