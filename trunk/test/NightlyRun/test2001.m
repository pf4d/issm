%Test Name: SquareSheetConstrainedGia2d
%GIA test, inspired on test101
md=triangle(model(),'../Exp/Square.exp',100000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');

%% indicate what you want to compute 
md.gia.cross_section_shape=1;    % for square-edged x-section 

%% define loading history 
md.timestepping.start_time=2400000; %2,400 kyr :: EVALUATION TIME
% to get rid of default final_time: make sure final_time>start_time
md.timestepping.final_time=2500000; %2,500 kyr
md.geometry.thickness=[...
	[md.geometry.thickness; 0.0],...
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
