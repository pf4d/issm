%Test Name: SquareShelfConstrainedMasstransp2dAdolc
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',1);
md.autodiff.isautodiff=true;
md.toolkits.DefaultAnalysis=issmgslsolver();
md.verbose=verbose('autodiff',true);
md=solve(md,'Masstransport');

%Fields and tolerances to track changes
field_names     ={'Thickness'};
field_tolerances={1e-13};
field_values={...
	(md.results.MasstransportSolution.Thickness),...
	};
