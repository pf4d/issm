import copy
from EnumDefinitions import *

def marshallcostfunctions(cost_functions):

	#copy list first
	data=copy.deepcopy(cost_functions)

	#convert to  Enums
	pos=[i for i,x in enumerate(cost_functions) if x==101];
	for i in pos: data[i]=SurfaceAbsVelMisfitEnum()        
	pos=[i for i,x in enumerate(cost_functions) if x==102];
	for i in pos: data[i]=SurfaceRelVelMisfitEnum()        
	pos=[i for i,x in enumerate(cost_functions) if x==103];
	for i in pos: data[i]=SurfaceLogVelMisfitEnum()        
	pos=[i for i,x in enumerate(cost_functions) if x==104];
	for i in pos: data[i]=SurfaceLogVxVyMisfitEnum()       
	pos=[i for i,x in enumerate(cost_functions) if x==105];
	for i in pos: data[i]=SurfaceAverageVelMisfitEnum()    
	pos=[i for i,x in enumerate(cost_functions) if x==201];
	for i in pos: data[i]=ThicknessAbsMisfitEnum()         
	pos=[i for i,x in enumerate(cost_functions) if x==501];
	for i in pos: data[i]=DragCoefficientAbsGradientEnum() 
	pos=[i for i,x in enumerate(cost_functions) if x==502];
	for i in pos: data[i]=RheologyBbarAbsGradientEnum()    
	pos=[i for i,x in enumerate(cost_functions) if x==503];
	for i in pos: data[i]=ThicknessAbsGradientEnum()       
	pos=[i for i,x in enumerate(cost_functions) if x==504];
	for i in pos: data[i]=ThicknessAlongGradientEnum()     
	pos=[i for i,x in enumerate(cost_functions) if x==505];
	for i in pos: data[i]=ThicknessAcrossGradientEnum()    

	return data
