import os
import numpy
from ContourToMesh import ContourToMesh

def SetIceSheetBC(md):
	"""
	SETICESHEETBC - Create the boundary conditions for stressbalance and thermal models for an IceSheet with no Ice Front

	   Usage:
	      md=SetIceSheetBC(md)

	   See also: SETICESHELFBC, SETMARINEICESHEETBC
	"""

	#node on Dirichlet
	pos=numpy.nonzero(md.mesh.vertexonboundary)
	md.stressbalance.spcvx=float('nan')*numpy.ones(md.mesh.numberofvertices)
	md.stressbalance.spcvy=float('nan')*numpy.ones(md.mesh.numberofvertices)
	md.stressbalance.spcvz=float('nan')*numpy.ones(md.mesh.numberofvertices)
	md.stressbalance.spcvx[pos]=0
	md.stressbalance.spcvy[pos]=0
	md.stressbalance.spcvz[pos]=0
	md.stressbalance.referential=float('nan')*numpy.ones((md.mesh.numberofvertices,6))
	md.stressbalance.loadingforce=0*numpy.ones((md.mesh.numberofvertices,3))

	#Dirichlet Values
	if isinstance(md.inversion.vx_obs,numpy.ndarray) and numpy.size(md.inversion.vx_obs,axis=0)==md.mesh.numberofvertices and isinstance(md.inversion.vy_obs,numpy.ndarray) and numpy.size(md.inversion.vy_obs,axis=0)==md.mesh.numberofvertices:
		print("      boundary conditions for stressbalance model: spc set as observed velocities")
		md.stressbalance.spcvx[pos]=md.inversion.vx_obs[pos]
		md.stressbalance.spcvy[pos]=md.inversion.vy_obs[pos]
	else:
		print("      boundary conditions for stressbalance model: spc set as zero")

	#No ice front -> do nothing

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
		if not isinstance(md.basalforcings.geothermalflux,numpy.ndarray) or not numpy.size(md.basalforcings.geothermalflux)==md.mesh.numberofvertices:
			md.basalforcings.geothermalflux=50.*10**-3*numpy.ones((md.mesh.numberofvertices,1))    #50 mW/m^2
	else:
		print("      no thermal boundary conditions created: no observed temperature found")

	return md

