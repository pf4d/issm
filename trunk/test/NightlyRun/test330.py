#Test Name:UnConfinedHydroDC
import numpy as np
import inspect
from model import *
from setmask import *
from triangle import triangle
from transient import transient
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from socket import gethostname
from generic import generic

md=triangle(model(),'../Exp/Strip.exp',10000.)
md=setmask(md,'','')
#reduced slab (20m long)
md.mesh.x=md.mesh.x/5.0e3
md.mesh.y=md.mesh.y/5.0e3
md=parameterize(md,'../Par/IceCube.py')
md.transient=transient.setallnullparameters(md.transient)
md.transient.ishydrology=True
md=setflowequation(md,'SSA','all')
md.cluster=generic('name',gethostname(),'np',1)
md.hydrology=hydrologydc()
md.hydrology=md.hydrology.initialize(md)

#Hydro Model Parameters
md.hydrology.isefficientlayer=0
md.hydrology.sedimentlimit_flag=0
md.hydrology.rel_tol=1.0e-6
md.hydrology.penalty_lock=0
md.hydrology.max_iter=200
md.hydrology.transfer_flag=0
md.hydrology.sediment_porosity=0.1
#Sediment
md.hydrology.sediment_thickness=10.0
md.hydrology.sediment_transmitivity=(1.0e-3*md.hydrology.sediment_thickness)*np.ones((md.mesh.numberofvertices))
#init
md.initialization.sediment_head=-5.0*np.ones((md.mesh.numberofvertices))
#BC
md.hydrology.spcsediment_head=np.nan*np.ones((md.mesh.numberofvertices))
md.hydrology.spcsediment_head[np.where(md.mesh.x==0)]=0.5

md.timestepping.time_step=5/md.constants.yts #5s steppin
md.settings.output_frequency=2
md.timestepping.final_time=300/md.constants.yts #500s run

md=solve(md,'Transient')

#fields to track, results can also be found in 
#Wang 2009 Fig 6b (jouranl of Hydrology)
field_names=['SedimentWaterHead1',
						 'SedimentWaterHead2']
field_tolerances=[1e-13, 
									1e-13]
field_values=[md.results.TransientSolution[10].SedimentHead,
							md.results.TransientSolution[30].SedimentHead]
