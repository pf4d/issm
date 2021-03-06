#Test Name: PigTranHO
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md=triangle(model(),'../Exp/Pig.exp',30000.)
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp')
md=parameterize(md,'../Par/Pig.py')
md.extrude(2,1.)
md=setflowequation(md,'HO','all')
md.cluster=generic('name',gethostname(),'np',3)
md=solve(md,'Transient')

# Fields and tolerances to track changes
field_names     =['Vx1','Vy1','Vz1','Vel1','Pressure1','Bed1','Surface1','Thickness1','Temperature1','BasalforcingsGroundediceMeltingRate1', \
				      'Vx2','Vy2','Vz2','Vel2','Pressure2','Bed2','Surface2','Thickness2','Temperature2','BasalforcingsGroundediceMeltingRate2']
field_tolerances=[1e-10,1e-10,1e-10,1e-10,1e-12,1e-11,2e-12,1e-11,1e-12,1e-09,\
						1e-11,1e-11,1e-09,1e-11,1e-11,1e-10,1e-11,1e-10,1e-11,2e-08]
field_values=[\
	md.results.TransientSolution[0].Vx,\
	md.results.TransientSolution[0].Vy,\
	md.results.TransientSolution[0].Vz,\
	md.results.TransientSolution[0].Vel,\
	md.results.TransientSolution[0].Pressure,\
	md.results.TransientSolution[0].Base,\
	md.results.TransientSolution[0].Surface,\
	md.results.TransientSolution[0].Thickness,\
	md.results.TransientSolution[0].Temperature,\
	md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate,\
	md.results.TransientSolution[1].Vx,\
	md.results.TransientSolution[1].Vy,\
	md.results.TransientSolution[1].Vz,\
	md.results.TransientSolution[1].Vel,\
	md.results.TransientSolution[1].Pressure,\
	md.results.TransientSolution[1].Base,\
	md.results.TransientSolution[1].Surface,\
	md.results.TransientSolution[1].Thickness,\
	md.results.TransientSolution[1].Temperature,\
	md.results.TransientSolution[1].BasalforcingsGroundediceMeltingRate,\
	]
