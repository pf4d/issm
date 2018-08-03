#Test Name: SquareSheetConstrainedGia2d
import numpy as np
import copy
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


#Define a model 
md=model()
md=triangle(md,'../Exp/Square.exp',100000.)
md=setmask(md,'','')
md=parameterize(md,'../Par/SquareSheetConstrained.py')

#Indicate what you want to compute 
md.gia.cross_section_shape=1    # for square-edged x-section 

#Define loading history (see test2001.m for the description)
md.timestepping.start_time=2400000 # 2,400 kyr
md.timestepping.final_time=2500000 # 2,500 kyr
md.geometry.thickness=np.vstack((np.hstack((md.geometry.thickness*0.0, 0.0)),
																 np.hstack((md.geometry.thickness/2.0, 0.1)),
																 np.hstack((md.geometry.thickness, 0.2)),
																 np.hstack((md.geometry.thickness, 1.0)),
																 np.hstack((md.geometry.thickness, md.timestepping.start_time)))).T

#Solve for GIA deflection 
md.cluster=generic('name',gethostname(),'np',3)
md=solve(md,'Gia')

#Fields and tolerances to track changes
field_names     =['GiaW','GiadWdt']
field_tolerances=[1e-13,1e-13]
field_values    =[md.results.GiaSolution.GiaW,
									md.results.GiaSolution.GiadWdt]

