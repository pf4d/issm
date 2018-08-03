import numpy
from SetMarineIceSheetBC import SetMarineIceSheetBC

#Ok, start defining model parameters here

print "      creating thickness"
ymin=numpy.min(md.mesh.y)
ymax=numpy.max(md.mesh.y)
md.geometry.thickness=500.*numpy.ones((md.mesh.numberofvertices))
md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness
md.geometry.surface=md.geometry.base+md.geometry.thickness

print "      creating drag"
md.friction.coefficient=200.*numpy.ones((md.mesh.numberofvertices))
md.friction.coefficient[numpy.nonzero(md.mask.groundedice_levelset<0.)[0]]=0.
md.friction.p=numpy.ones((md.mesh.numberofelements))
md.friction.q=numpy.ones((md.mesh.numberofelements))

print "      creating initial values"
md.initialization.temperature=(273.-20.)*numpy.ones((md.mesh.numberofvertices))
md.initialization.vx=numpy.zeros((md.mesh.numberofvertices))
md.initialization.vy=numpy.zeros((md.mesh.numberofvertices))
md.initialization.vz=numpy.zeros((md.mesh.numberofvertices))
md.initialization.vel=numpy.zeros((md.mesh.numberofvertices))
md.initialization.pressure=numpy.zeros((md.mesh.numberofvertices))

print "      creating flow law parameter"
md.materials.rheology_B=1.7687*10**8*numpy.ones((md.mesh.numberofvertices))
md.materials.rheology_n=3.*numpy.ones((md.mesh.numberofelements))

print "      creating surface mass balance"
md.smb.mass_balance=0.2*numpy.ones((md.mesh.numberofvertices))    #0m/a
md.basalforcings.floatingice_melting_rate=0.*numpy.ones((md.mesh.numberofvertices))    #0m/a
md.basalforcings.groundedice_melting_rate=0.*numpy.ones((md.mesh.numberofvertices))    #0m/a

print "      boundary conditions"
md=SetMarineIceSheetBC(md,'../Exp/SquareFrontEISMINT.exp')

#Evolution of the ice shelf
pos=numpy.nonzero(md.mesh.y==200000.)    #nodes on the upper boundary condition
md.balancethickness.spcthickness=float('NaN')*numpy.ones((md.mesh.numberofvertices))
md.balancethickness.spcthickness[pos]=500.
md.masstransport.spcthickness=float('NaN')*numpy.ones((md.mesh.numberofvertices))
md.masstransport.spcthickness[pos]=500.
md.masstransport.stabilization=0    #Better result with no artificial diffusivity
md.thermal.stabilization=0
md.timestepping.final_time=500.
md.timestepping.time_step=1
