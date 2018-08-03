import os
import numpy
from ContourToMesh import ContourToMesh
import MatlabFuncs as m

def SetMarineIceSheetBC(md,icefrontfile=''):
	"""
	SETICEMARINESHEETBC - Create the boundary conditions for stressbalance and thermal models for a  Marine Ice Sheet with Ice Front

	   Neumann BC are used on the ice front (an ARGUS contour around the ice front
	   can be given in input, or it will be deduced as onfloatingice & onboundary)
	   Dirichlet BC are used elsewhere for stressbalance

	   Usage:
	      md=SetMarineIceSheetBC(md,icefrontfile)
	      md=SetMarineIceSheetBC(md)

	   Example:
	      md=SetMarineIceSheetBC(md,'Front.exp')
	      md=SetMarineIceSheetBC(md)

	   See also: SETICESHELFBC, SETMARINEICESHEETBC
	"""

	#node on Dirichlet (boundary and ~icefront)
	if icefrontfile:
		#User provided Front.exp, use it
		if not os.path.exists(icefrontfile):
			raise IOError("SetMarineIceSheetBC error message: ice front file '%s' not found." % icefrontfile)
		[incontour,dum]=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,icefrontfile,'node',2)
		vertexonicefront=numpy.logical_and(md.mesh.vertexonboundary,incontour.reshape(-1))
	else:
		#Guess where the ice front is
		vertexonfloatingice=numpy.zeros((md.mesh.numberofvertices,1))
		pos=numpy.nonzero(numpy.sum(md.mask.groundedice_levelset[md.mesh.elements-1]<0.,axis=1) >0.)[0]
		vertexonfloatingice[md.mesh.elements[pos].astype(int)-1]=1.
		vertexonicefront=numpy.logical_and(numpy.reshape(md.mesh.vertexonboundary,(-1,1)),vertexonfloatingice>0.)

#	pos=find(md.mesh.vertexonboundary & ~vertexonicefront);
	pos=numpy.nonzero(numpy.logical_and(md.mesh.vertexonboundary,numpy.logical_not(vertexonicefront)))[0]
	if not numpy.size(pos):
		print("SetMarineIceSheetBC warning: ice front all around the glacier, no dirichlet found. Dirichlet must be added manually.")

	md.stressbalance.spcvx=float('nan')*numpy.ones(md.mesh.numberofvertices)
	md.stressbalance.spcvy=float('nan')*numpy.ones(md.mesh.numberofvertices)
	md.stressbalance.spcvz=float('nan')*numpy.ones(md.mesh.numberofvertices)
	md.stressbalance.referential=float('nan')*numpy.ones((md.mesh.numberofvertices,6))
	md.stressbalance.loadingforce=0*numpy.ones((md.mesh.numberofvertices,3))

	#Position of ice front
	pos=numpy.nonzero(vertexonicefront)[0]
	md.mask.ice_levelset[pos]=0

	#First find segments that are not completely on the front
	if m.strcmp(md.mesh.elementtype(),'Penta'):
		numbernodesfront=4
	elif m.strcmp(md.mesh.elementtype(),'Tria'):
		numbernodesfront=2
	else:
			raise Exception("Mesh type not supported")
	if any(md.mask.ice_levelset<=0):
		values=md.mask.ice_levelset[md.mesh.segments[:,0:-1]-1]
		segmentsfront=1-values
		numpy.sum(segmentsfront,axis=1)!=numbernodesfront
		segments=numpy.nonzero(numpy.sum(segmentsfront,axis=1)!=numbernodesfront)[0]
		#Find all nodes for these segments and spc them
		pos=md.mesh.segments[segments,0:-1]-1
	else:
		pos=numpy.nonzero(md.mesh.vertexonboundary)[0]
	md.stressbalance.spcvx[pos]=0
	md.stressbalance.spcvy[pos]=0
	md.stressbalance.spcvz[pos]=0

	#Dirichlet Values
	if isinstance(md.inversion.vx_obs,numpy.ndarray) and numpy.size(md.inversion.vx_obs,axis=0)==md.mesh.numberofvertices and isinstance(md.inversion.vy_obs,numpy.ndarray) and numpy.size(md.inversion.vy_obs,axis=0)==md.mesh.numberofvertices:
		print("      boundary conditions for stressbalance model: spc set as observed velocities")
		md.stressbalance.spcvx[pos]=md.inversion.vx_obs[pos]
		md.stressbalance.spcvy[pos]=md.inversion.vy_obs[pos]
	else:
		print("      boundary conditions for stressbalance model: spc set as zero")

	md.hydrology.spcwatercolumn=numpy.zeros((md.mesh.numberofvertices,2))
	pos=numpy.nonzero(md.mesh.vertexonboundary)[0]
	md.hydrology.spcwatercolumn[pos,0]=1

	#Create zeros basalforcings and smb
	md.smb.initialize(md)
	md.basalforcings.initialize(md)

	#Deal with other boundary conditions
	if numpy.all(numpy.isnan(md.balancethickness.thickening_rate)):
		md.balancethickness.thickening_rate=numpy.zeros((md.mesh.numberofvertices,1))
		print("      no balancethickness.thickening_rate specified: values set as zero")

	md.masstransport.spcthickness=float('nan')*numpy.ones((md.mesh.numberofvertices,1))
	md.balancethickness.spcthickness=float('nan')*numpy.ones((md.mesh.numberofvertices,1))
	md.damage.spcdamage=float('nan')*numpy.ones((md.mesh.numberofvertices,1))

	if isinstance(md.initialization.temperature,numpy.ndarray) and numpy.size(md.initialization.temperature,axis=0)==md.mesh.numberofvertices:
		md.thermal.spctemperature=float('nan')*numpy.ones((md.mesh.numberofvertices,1))
		if hasattr(md.mesh,'vertexonsurface'):
			pos=numpy.nonzero(md.mesh.vertexonsurface)[0]
			md.thermal.spctemperature[pos]=md.initialization.temperature[pos]    #impose observed temperature on surface
		if not isinstance(md.basalforcings.geothermalflux,numpy.ndarray) or not numpy.size(md.basalforcings.geothermalflux,axis=0)==md.mesh.numberofvertices:
			md.basalforcings.geothermalflux=numpy.zeros((md.mesh.numberofvertices,1))
			md.basalforcings.geothermalflux[numpy.nonzero(md.mask.groundedice_levelset>0.)]=50.*10.**-3    #50mW/m2
	else:
		print("      no thermal boundary conditions created: no observed temperature found")

	return md

