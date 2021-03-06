import numpy
import os
import model
import glob
def exportVTK(filename,model,*args):
	'''
	vtk export
	function exportVTK(filename,model)
	creates a directory with the vtk files for displays in paraview
	(only work for triangle and wedges based on their number of nodes)
	
	Give only the results for nw but could be extended to geometry, mask... 
	
	input: filename   destination 
	(string)
	------------------------------------------------------------------
model      this is md 
	------------------------------------------------------------------
	By default only the results are exported, you can add whichever
	field you need as a string:
	add 'geometry' to export md.geometry

	Basile de Fleurian:
	'''
	Dir=os.path.basename(filename)
	Path=filename[:-len(Dir)]

	if os.path.exists(filename):
		print(('File {} allready exist'.format(filename)))
		newname=input('Give a new name or "delete" to replace: ')
		if newname=='delete':
			filelist = glob.glob(filename+'/*')
			for oldfile in filelist:
				os.remove(oldfile)
		else:
			print(('New file name is {}'.format(newname)))
			filename=newname
			os.mkdir(filename)
	else:
		os.mkdir(filename)

	#get the element related variables
	if 'z' in dict.keys(model.mesh.__dict__):
		points=numpy.column_stack((model.mesh.x,model.mesh.y,model.mesh.z))
		dim=3
	else:
		points=numpy.column_stack((model.mesh.x,model.mesh.y,numpy.zeros(numpy.shape(model.mesh.x))))
		dim=2

	num_of_points=numpy.size(model.mesh.x)
	num_of_elt=numpy.shape(model.mesh.elements)[0]
	point_per_elt=numpy.shape(model.mesh.elements)[1]
		
	#Select the type of element function of the number of nodes per elements
	if point_per_elt==3:
		celltype=5 #triangles
	elif point_per_elt==6:
		celltype=13 #wedges
	else:
		error('Your Element definition is not taken into account \n')

	#this is the result structure
	res_struct=model.results
	if (len(res_struct.__dict__)>0):
		#Getting all the solutions of the model
		solnames=(dict.keys(res_struct.__dict__))
		num_of_sols=len(solnames)
		num_of_timesteps=1
		out_freq=model.settings.output_frequency
		#%building solutionstructure 
		for solution in solnames:
			#looking for multiple time steps
			if (numpy.size(res_struct.__dict__[solution])>num_of_timesteps):
				num_of_timesteps=numpy.size(res_struct.__dict__[solution])
				num_of_timesteps=int(num_of_timesteps/out_freq)+1
	else:
		num_of_timesteps=1

	for step in range(0,num_of_timesteps):
		timestep=step
		fid=open((filename +'/Timestep.vtk'+str(timestep)+'.vtk'),'w+')
		fid.write('# vtk DataFile Version 2.0 \n')
		fid.write('Data for run %s \n' % model.miscellaneous.name)
		fid.write('ASCII \n')
		fid.write('DATASET UNSTRUCTURED_GRID \n')
		fid.write('POINTS %d float\n' % num_of_points)
		if(dim==3):
			for point in points:
				fid.write('%f %f %f \n'%(point[0], point[1], point[2]))
		elif(dim==2):
			for point in points:
				fid.write('%f %f %f \n'%(point[0], point[1], point[2]))
			
		fid.write('CELLS %d %d\n' %(num_of_elt, num_of_elt*(point_per_elt+1)))
		
		if point_per_elt==3:
			for elt in range(0, num_of_elt):
				fid.write('3 %d %d %d\n' %(model.mesh.elements[elt,0]-1,model.mesh.elements[elt,1]-1,model.mesh.elements[elt,2]-1))
		elif point_per_elt==6:
			for elt in range(0, num_of_elt):
				fid.write('6 %d %d %d %d %d %d\n' %(model.mesh.elements[elt,0]-1,model.mesh.elements[elt,1]-1,model.mesh.elements[elt,2]-1,model.mesh.elements[elt,3]-1,model.mesh.elements[elt,4]-1,model.mesh.elements[elt,5]-1))
		else:
			print('Number of nodes per element not supported')

		fid.write('CELL_TYPES %d\n' %num_of_elt)
		for elt in range(0, num_of_elt):
			fid.write('%d\n' %celltype)

		fid.write('POINT_DATA %s \n' %str(num_of_points))
	
		#loop over the different solution structures
		if 'solnames' in locals():
			for sol in solnames:
				#dealing with results on different timesteps
				if(numpy.size(res_struct.__dict__[sol])>timestep):
					timestep = step
				else:
					timestep = numpy.size(res_struct.__dict__[sol])
				
				#getting the  fields in the solution
				if(numpy.size(res_struct.__dict__[sol])>1):
					fieldnames=dict.keys(res_struct.__dict__[sol].__getitem__(timestep*out_freq-1).__dict__)
				else:
					fieldnames=dict.keys(res_struct.__dict__[sol].__dict__)
				#check which field is a real result and print
				for field in fieldnames:
					if(numpy.size(res_struct.__dict__[sol])>1):
						fieldstruct=res_struct.__dict__[sol].__getitem__(timestep*out_freq-1).__dict__[field]
					else:
						fieldstruct=res_struct.__dict__[sol].__dict__[field]

					if ((numpy.size(fieldstruct))==num_of_points):
						fid.write('SCALARS %s float 1 \n' % field)
						fid.write('LOOKUP_TABLE default\n')
						for node in range(0,num_of_points):
							#paraview does not like NaN, replacing
							if numpy.isnan(fieldstruct[node]):
								fid.write('%e\n' % -9999.9999)
							#also checking for verry small value that mess up
							elif (abs(fieldstruct[node])<1.0e-20):
								fid.write('%e\n' % 0.0)
							else:
								fid.write('%e\n' % fieldstruct[node])
					
		#loop on arguments, if something other than result is asked, do
		#it now
		for other in args:
			other_struct=model.__dict__[other]
			othernames=(dict.keys(other_struct.__dict__))
			for field in othernames:
				if ((numpy.size(other_struct.__dict__[field]))==num_of_points):
					fid.write('SCALARS %s float 1 \n' % field)
					fid.write('LOOKUP_TABLE default\n')
					for node in range(0,num_of_points):
						#paraview does not like NaN, replacing
						if numpy.isnan(other_struct.__dict__[field][node]):
							fid.write('%e\n' % -9999.9999)
						#also checking for verry small value that mess up
						elif (abs(other_struct.__dict__[field][node])<1.0e-20):
							fid.write('%e\n' % 0.0)
						else:
							fid.write('%e\n' % other_struct.__dict__[field][node])
	fid.close();
