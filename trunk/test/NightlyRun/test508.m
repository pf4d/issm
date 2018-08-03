%Test Name: PigSteaSSA3d
md=triangle(model(),'../Exp/Pig.exp',35000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md=extrude(md,3,1.1);
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md.timestepping.time_step=0.;
md=solve(md,'Steadystate');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vz','Vel','Pressure','Temperature','BasalforcingsGroundediceMeltingRate'};
field_tolerances={5e-08,3e-08,5e-08,3e-08,1e-09,2e-07,8e-07};
field_values={...
	(md.results.SteadystateSolution.Vx),...
	(md.results.SteadystateSolution.Vy),...
	(md.results.SteadystateSolution.Vz),...
	(md.results.SteadystateSolution.Vel),...
	(md.results.SteadystateSolution.Pressure),...
	(md.results.SteadystateSolution.Temperature),...
	(md.results.SteadystateSolution.BasalforcingsGroundediceMeltingRate),...
	};