%Test Name: GiaIvinsBenchmarksAB2dC3
% Benchmark experiments (Figure A2c Ivins and James, 1999, Geophys. J. Int.) 
md=triangle(model(),'../Exp/RoundFrontEISMINT.exp',200000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/GiaIvinsBenchmarksCD.par');

%% indicate what you want to compute 
md.gia.cross_section_shape=1;    % for square-edged x-section 

%% define loading history 
md.timestepping.start_time=2400000; % for t \approx \infty 
md.timestepping.final_time=2500000; % 2,500 kyr
md.geometry.thickness=[...
	[md.geometry.thickness*0.0; 0.0],...
	[md.geometry.thickness/2.0; 0.1],...
	[md.geometry.thickness; 0.2],...
	[md.geometry.thickness; md.timestepping.start_time],...
	];

%% solve for GIA deflection 
md.cluster=generic('name',oshostname(),'np',3);
md.verbose=verbose('1111111');
md=solve(md,'Gia');

%Fields and tolerances to track changes
field_names     ={'GiaW','GiadWdt'};
field_tolerances={1e-13,1e-13};
field_values={...
	(md.results.GiaSolution.GiaW),...
	(md.results.GiaSolution.GiadWdt),...
	};
