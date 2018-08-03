#Test Name: SquareSheetShelfStressHO
import numpy as np
from model import *
from socket import gethostname

from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md=triangle(model(),'../Exp/Square.exp',180000.)
md=setmask(md,'../Exp/SquareShelf.exp','')
md=parameterize(md,'../Par/SquareSheetShelf.py')
md.extrude(5,1.)
md=setflowequation(md,'HO','all')
md.cluster=generic('name',gethostname(),'np',3)
md=solve(md,'Stressbalance')

#Fields and tolerances to track changes
field_names     =['Vx','Vy','Vz','Vel','Pressure']
field_tolerances=[2e-09,2e-09,2e-09,2e-09,2e-09]
field_values=[\
	md.results.StressbalanceSolution.Vx,\
	md.results.StressbalanceSolution.Vy,\
	md.results.StressbalanceSolution.Vz,\
	md.results.StressbalanceSolution.Vel,\
	md.results.StressbalanceSolution.Pressure,\
	]
