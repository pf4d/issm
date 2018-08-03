%Test Name: PigSteaHO
md=triangle(model(),'../Exp/Pig.exp',30000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md=extrude(md,3,1.);
md=setflowequation(md,'HO','all');
md.cluster=generic('name',oshostname(),'np',3);
md.timestepping.time_step=0.;
md.thermal.penalty_threshold=7;
md=solve(md,'Steadystate');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vz','Vel','Pressure','Temperature','BasalforcingsGroundediceMeltingRate'};
field_tolerances={1e-09,1e-09,1e-09,1e-09,1e-09,1e-09,1e-06
};
field_values={...
	(md.results.SteadystateSolution.Vx),...
	(md.results.SteadystateSolution.Vy),...
	(md.results.SteadystateSolution.Vz),...
	(md.results.SteadystateSolution.Vel),...
	(md.results.SteadystateSolution.Pressure),...
	(md.results.SteadystateSolution.Temperature),...
	(md.results.SteadystateSolution.BasalforcingsGroundediceMeltingRate),...
	};