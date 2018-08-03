#Test Name: SquareSheetConstrainedHydrologyDC
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from transient import *
from setflowequation import *
from solve import *


from generic import generic

md=triangle(model(),'../Exp/Square.exp',100000.)
md=setmask(md,'','')
md=parameterize(md,'../Par/IceCube.py')

md.transient=transient.setallnullparameters(md.transient)
md.transient.ishydrology=True

md=setflowequation(md,'SSA','all')
md.cluster=generic('name',gethostname(),'np',1)
md.hydrology=hydrologydc()
md.hydrology=md.hydrology.initialize(md)

md.hydrology.isefficientlayer=0
md.hydrology.sedimentlimit_flag=1
md.hydrology.sedimentlimit=8000.0
md.initialization.sediment_head=np.zeros((md.mesh.numberofvertices))
md.hydrology.spcsediment_head=np.nan*np.ones((md.mesh.numberofvertices))
pos=np.nonzero(md.mesh.y==0.)[0]
md.hydrology.spcsediment_head[pos]=0.0
md.basalforcings.groundedice_melting_rate = 2.0*np.ones((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate = 0.0*np.ones((md.mesh.numberofvertices))
md.hydrology.sediment_transmitivity= 3.0*np.ones((md.mesh.numberofvertices))
md.timestepping.time_step=0
md.timestepping.final_time=1.0
md=solve(md,'Hydrology')

#Fields and tolerances to track changes
#you can also compare with an analitic solution, but it is exact
#only if no limits are applied
#analitic=(md.mesh.y**2-2*md.mesh.y*1.0e6)*(-2.0/(2*md.constants.yts*md.hydrology.sediment_transmitivity))
field_names     =['SedimentWaterHead','SedimentHeadResidual']
field_tolerances=[1e-13, 3e-10]
field_values=[md.results.HydrologySolution.SedimentHead,md.results.HydrologySolution.SedimentHeadResidual]
