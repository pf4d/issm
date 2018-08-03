#Test Name: SquareShelfStressFS

from model import *
from socket import gethostname
import numpy as np
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from generic import generic

md=triangle(model(),'../Exp/Square.exp',180000)
md=setmask(md,'all','')
md=parameterize(md,'../Par/SquareShelf.py')
md.extrude(3,2.)
md=setflowequation(md,'FS','all')
md.cluster=generic('name',gethostname(),'np',3)
md.stressbalance.shelf_dampening=1
md.timestepping.time_step=0
md1=solve(md,'Stressbalance')
md.stressbalance.shelf_dampening=0
md=solve(md,'Stressbalance')


# Fields and tolerances to track changes

field_names     =['Vx','Vy','Vz','Vel','Pressure']
field_tolerances=[1e-08,1e-08,2e-06,1e-08,1e-08]
field_values=[\
	md.results.StressbalanceSolution.Vx,\
	md.results.StressbalanceSolution.Vy,\
	md.results.StressbalanceSolution.Vz,\
	md.results.StressbalanceSolution.Vel,\
	md.results.StressbalanceSolution.Pressure,\
	]
