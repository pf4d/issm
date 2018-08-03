%Test Name: SquareSheetConstrainedEnthalpyTran
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');
md=extrude(md,3,1.);
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md.initialization.waterfraction=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);
md.transient.isstressbalance=0;
md.transient.ismasstransport=0;
md.transient.issmb=1;
md.transient.isthermal=1;
md.transient.isgroundingline=0;
md.thermal.isenthalpy=1;
md.thermal.isdynamicbasalspc=1;
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'Enthalpy1','Waterfraction1','Temperature1',...
	'Enthalpy2','Waterfraction2','Temperature2',...
	'Enthalpy3','Waterfraction3','Temperature3'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.TransientSolution(1).Enthalpy),...
	(md.results.TransientSolution(1).Waterfraction),...
	(md.results.TransientSolution(1).Temperature),...
	(md.results.TransientSolution(2).Enthalpy),...
	(md.results.TransientSolution(2).Waterfraction),...
	(md.results.TransientSolution(2).Temperature),...
	(md.results.TransientSolution(3).Enthalpy),...
	(md.results.TransientSolution(3).Waterfraction),...
	(md.results.TransientSolution(3).Temperature),...
	};