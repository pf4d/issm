from netCDF4 import Dataset, stringtochar
import numpy as np
import time
import collections
from mesh2d import *
from mesh3dprisms import *
from results import *
from os import path, remove

def export_netCDF(md,filename):
	#Now going on Real treatment
	if path.exists(filename):
		print ('File {} allready exist'.format(filename))
		newname=raw_input('Give a new name or "delete" to replace: ')
		if newname=='delete':
			remove(filename)
		else:
			print ('New file name is {}'.format(newname))
			filename=newname
			
	NCData=Dataset(filename, 'w', format='NETCDF4')
	NCData.description = 'Results for run' + md.miscellaneous.name
	NCData.history = 'Created ' + time.ctime(time.time())

	#define netCDF dimensions
	try:
		StepNum=np.shape(dict.values(md.results.__dict__))[1]
	except IndexError:
		StepNum=1
	Dimension1=NCData.createDimension('DimNum1',StepNum)#time is first
	DimDict={len(Dimension1):'DimNum1'}
	dimindex=1

	dimlist=[2,md.mesh.numberofelements,md.mesh.numberofvertices,np.shape(md.mesh.elements)[1]]
	for i in range(0,4):
		if dimlist[i] not in DimDict.keys():
			dimindex+=1
			NewDim=NCData.createDimension('DimNum'+str(dimindex),dimlist[i])
			DimDict[len(NewDim)]='DimNum'+str(dimindex)

	typelist=[bool,str,unicode,int,float,complex,
						collections.OrderedDict,
						np.int64,np.ndarray,np.float64]
	groups=dict.keys(md.__dict__)
	#get all model classes and create respective groups
	for group in groups:
		NCgroup=NCData.createGroup(str(group))
		#In each group gather the fields of the class
		fields=dict.keys(md.__dict__[group].__dict__)

		#looping on fields
		for field in fields:
			#Special treatment for list fields
			if type(md.__dict__[group].__dict__[field])==list:
				StdList=False
				if len(md.__dict__[group].__dict__[field])==0:
					StdList=True
				else:
					StdList=type(md.__dict__[group].__dict__[field][0]) in typelist
				NCgroup.__setattr__('classtype', md.__dict__[group].__class__.__name__)
				if StdList: #this is a standard or empty list just proceed
					Var=md.__dict__[group].__dict__[field]
					DimDict=CreateVar(NCData,Var,field,NCgroup,DimDict)
				else: #this is a list of fields, specific treatment needed
					Listsize=len(md.__dict__[group].__dict__[field])
					Subgroup=NCgroup.createGroup(str(field))
					Subgroup.__setattr__('classtype',md.__dict__[group].__dict__[field].__class__.__name__)
					for listindex in range(0,Listsize):
						try:
							Listgroup=Subgroup.createGroup(str(md.__dict__[group].__dict__[field].__getitem__(listindex).__dict__['name']))
						except KeyError:
							for naming in ['step']:
								Listgroup=Subgroup.createGroup(str(md.__dict__[group].__dict__[field].__getitem__(listindex).__dict__[naming]))
						except AttributeError:
							Listgroup=Subgroup.createGroup(str(md.__dict__[group].__dict__[field].__class__.__name__)+str(listindex))
						Listgroup.__setattr__('classtype',md.__dict__[group].__dict__[field].__getitem__(listindex).__class__.__name__)
						try:
							subfields=dict.keys(md.__dict__[group].__dict__[field].__getitem__(listindex).__dict__)
						except AttributeError:
							subfields=dict.keys(md.__dict__[group].__dict__[field].__getitem__(listindex))
						for subfield in subfields:
							if subfield!='outlog':
								try:
									Var=md.__dict__[group].__dict__[field].__getitem__(listindex).__dict__[subfield]
								except AttributeError:
									Var=md.__dict__[group].__dict__[field].__getitem__(listindex)[subfield]
								DimDict=CreateVar(NCData,Var,subfield,Listgroup,DimDict,md.__dict__[group],field,listindex)

			#No subgroup, we directly treat the variable
			elif type(md.__dict__[group].__dict__[field]) in typelist or field=='bamg':
				NCgroup.__setattr__('classtype', md.__dict__[group].__class__.__name__)
				Var=md.__dict__[group].__dict__[field]
				DimDict=CreateVar(NCData,Var,field,NCgroup,DimDict)
			elif md.__dict__[group].__dict__[field] is None:
				print( 'field md.{}.{} is None'.format(group,field))
				#do nothing
			else:
				NCgroup.__setattr__('classtype', str(group))
				Subgroup=NCgroup.createGroup(str(field))
				Subgroup.__setattr__('classtype',md.__dict__[group].__class__.__name__)
				subfields=dict.keys(md.__dict__[group].__dict__[field].__dict__)
                                
				for subfield in subfields:
					if str(subfield)!='outlog':
						Var=md.__dict__[group].__dict__[field].__dict__[subfield]
						DimDict=CreateVar(NCData,Var,subfield,Subgroup,DimDict)
				
	NCData.close()

#============================================================================
#Define the variables
def CreateVar(NCData,var,field,Group,DimDict,*step_args):
	#grab type
	try:
		val_type=str(var.dtype)
	except AttributeError:
		val_type=type(var)
		#grab dimension
	try:
		val_shape=dict.keys(var)
	except TypeError:
		val_shape=np.shape(var)

	TypeDict = {float:'f8',
							'float64':'f8',
							np.float64:'f8',
							int:'i8',
							'int64':'i8',
							np.int64:'i8',
							str:str,
							dict:str}
		
	val_dim=np.shape(val_shape)[0]
	#Now define and fill up variable
	#treating scalar string or bool as atribute
	if val_type==str or val_type==unicode or val_type==bool:
		Group.__setattr__(str(field).swapcase(), str(var))
	#treating list as string table
	elif val_type==list:
		dimensions,DimDict=GetDim(NCData,var,val_shape,DimDict,val_dim)
		#try to get the type from the first element
		try:
			nctype=TypeDict[type(var[0])]
		except IndexError:
			nctype=str #most probably an empty list take str for that
		ncvar = Group.createVariable(str(field),nctype,dimensions,zlib=True)
		if val_shape==0:
			ncvar= []
		else:			
			for elt in range(0,val_shape[0]):
				ncvar[elt] = var[elt]
	#treating bool tables as string tables
	elif val_type=='bool':
		dimensions,DimDict=GetDim(NCData,var,val_shape,DimDict,val_dim)
		ncvar = Group.createVariable(str(field),str,dimensions,zlib=True)
		for elt in range(0,val_shape[0]):
			ncvar[elt] = str(var[elt])
	#treating dictionaries as tables of strings
	elif val_type==collections.OrderedDict or val_type==dict:
		dimensions,DimDict=GetDim(NCData,var,val_shape,DimDict,val_dim)
		ncvar = Group.createVariable(str(field),str,dimensions,zlib=True)
		for elt in range(0,val_dim):
			ncvar[elt,0]=dict.keys(var)[elt]
			ncvar[elt,1]=str(dict.values(var)[elt]) #converting to str to avoid potential problems
	#Now dealing with numeric variables
	else:
		dimensions,DimDict=GetDim(NCData,var,val_shape,DimDict,val_dim)
		ncvar = Group.createVariable(str(field),TypeDict[val_type],dimensions,zlib=True)
		try:
			nan_val=np.isnan(var)
			if nan_val.all():
				ncvar [:] = 'NaN'
			else:
				ncvar[:] = var
		except TypeError: #type does not accept nan, get vallue of the variable
			ncvar[:] = var
	return DimDict

#============================================================================
#retriev the dimension tuple from a dictionnary
def GetDim(NCData,var,shape,DimDict,i):
	output=[]
	#grab dimension
	for dim in range(0,i): #loop on the dimensions
		if type(shape[0])==int:
			try:
				output=output+[str(DimDict[shape[dim]])] #test if the dimension allready exist
			except KeyError: #if not create it
				if (shape[dim])>0:
					index=len(DimDict)+1
					NewDim=NCData.createDimension('DimNum'+str(index),(shape[dim]))
					DimDict[len(NewDim)]='DimNum'+str(index)
					output=output+[str(DimDict[shape[dim]])]
		elif type(shape[0])==str or type(shape[0])==unicode:#dealling with a dictionnary
			try:
				#dimension5 is 2 to treat with dict
				output=[str(DimDict[np.shape(shape)[0]])]+[DimDict[2]]
			except KeyError:
				index=len(DimDict)+1
				NewDim=NCData.createDimension('DimNum'+str(index),np.shape(shape)[0])
				DimDict[len(NewDim)]='DimNum'+str(index)
				output=[str(DimDict[np.shape(dict.keys(var))[0]])]+[DimDict[2]]
			break
	return tuple(output), DimDict
