from math import isnan
import numpy as np

def processdata(md,data,options):
    """
    PROCESSDATA - process data to be plotted
    
    	datatype = 1 -> elements
    	datatype = 2 -> nodes
    	datatype = 3 -> node quivers
    	datatype = 4 -> patch
    
    	Usage:
    		data,datatype=processdata(md,data,options);
    
    	See also: PLOTMODEL, PROCESSMESH
    """
    
    #check format
    if (len(data)==0 or (len(data)==1 and not isinstance(data,dict) and isnan(data).all())):
        raise ValueError("processdata error message: 'data' provided is empty")
    
    #needed later on
    if 'numberofvertices2d' in dir(md.mesh):
    	numberofvertices2d=md.mesh.numberofvertices2d
    	numberofelements2d=md.mesh.numberofelements2d
    else:
    	numberofvertices2d=np.nan
    	numberofelements2d=np.nan
    
    procdata=np.copy(data)
    
    #process patch
    
    #initialize datatype
    datatype=0
    
    #get datasize
    if np.ndim(procdata)==1:
    	datasize=np.array([len(procdata),1])
    else:
    	datasize=np.shape(procdata)
        if len(datasize)>2:
            raise ValueError('data passed to plotmodel has more than 2 dimensions; check that column vectors are rank-1')
    
    #process NaN's if any
    nanfill=options.getfieldvalue('nan',-9999)
    if np.any(np.isnan(procdata)):
    	lb=np.min(data[~np.isnan(data)])
    	ub=np.max(data[~np.isnan(data)])
    	if lb==ub:
    	    lb=lb-0.5
    	    ub=ub+0.5
    	    nanfill=lb-1
    	procdata[np.isnan(procdata)]=nanfill
    	options.addfielddefault('clim',[lb,ub])
    	options.addfielddefault('cmap_set_under','1')
    	print(("WARNING: nan's treated as", nanfill, "by default.  Change using pairoption 'nan',nan_fill_value in plotmodel call"))
    
    #quiver plot
    if datasize[1]>1 and datasize[0]!= md.mesh.numberofvertices+1:
        if datasize[0]==md.mesh.numberofvertices and datasize[1]==2:
            datatype=3
        else:
            raise ValueError('plotmodel error message: data should have two columns of length md.mesh.numberofvertices for a quiver plot')
    
    #non-patch processing 
    
    #element data
    if datasize[0]==md.mesh.numberofelements and datasize[1]==1:
    	
    	#initialize datatype if non patch
    	if datatype!=4 and datatype!=5:
    	    datatype=1
    
    	#mask?
    
    	#log?
    
    #node data
    if datasize[0]==md.mesh.numberofvertices and datasize[1]==1:
    	datatype=2
    
    #spc time series? 
    if datasize[0]==md.mesh.numberofvertices+1:
    	datatype=2
        spccol=options.getfieldvalue('spccol',0)
        print('multiple-column spc field; specify column to plot using option "spccol"')
        print(('column ', spccol, ' plotted for time: ', procdata[-1,spccol]))
        procdata=procdata[0:-1,spccol]
    
    	#mask?
    
    	#log?
    
    #layer projection?
    
    #control arrow density if quiver plot
    
    #convert rank-2 array to rank-1
    if np.ndim(procdata)==2 and np.shape(procdata)[1]==1:
    	procdata=procdata.reshape(-1,)
    
    #if datatype is still zero, error out
    if datatype==0:
    	raise ValueError("processdata error: data provided not recognized or not supported")
    else:
    	return procdata, datatype
