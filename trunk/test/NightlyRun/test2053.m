%Test Name: GiaIvinsBenchmarksAB2dA3
% Benchmark experiments (Figure A2a Ivins and James, 1999, Geophys. J. Int.) 
md=triangle(model(),'../Exp/RoundFrontEISMINT.exp',200000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/GiaIvinsBenchmarksAB.par');

%% indicate what you want to compute 
md.gia.cross_section_shape=1;    % for square-edged x-section 

%% define loading history 
md.timestepping.start_time=2010100; % after 10 kyr of deglaciation 
md.timestepping.final_time=2500000; % 2,500 kyr
md.geometry.thickness=[...
	[md.geometry.thickness*0.0; 0.0],...
	[md.geometry.thickness; 1000],...
	[md.geometry.thickness; 2000000],...
	[md.geometry.thickness*0.0; 2000100],...
	[md.geometry.thickness*0.0; md.timestepping.start_time],...
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
