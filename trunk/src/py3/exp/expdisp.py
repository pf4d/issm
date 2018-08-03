from expread import expread
import numpy as np

def expdisp(domainoutline,ax,linestyle='--k',linewidth=1,unitmultiplier=1.):
    '''
    plot the contents of a domain outline file

    This routine reads in a domain outline file and plots all of the x,y contours

    'ax' is a handle to the current plot axes, onto which the contours are to be drawn

    Usage:
        expdisp(domainoutline,ax)

    Example:
        expdisp('domain.exp',plt.gca(),linestyle='--k',linewidth=2,unitmultiplier=1.e3)
    '''

    domain=expread(domainoutline)

    for i in range(len(domain)):
        if domain[i]['nods']==1:
            ax.plot(domain[i]['x']*unitmultiplier,domain[i]['y']*unitmultiplier,'o',mec='k',mfc='r',ms=10)
        else:
            x=domain[i]['x'].tolist() # since expread returns a string representation of the arrays
            y=domain[i]['y'].tolist()
            ax.plot(x*unitmultiplier,y*unitmultiplier,linestyle,linewidth=linewidth)
