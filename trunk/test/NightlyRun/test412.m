%Test Name: SquareSheetShelfDiadSSA3dDakota
md=triangle(model(),'../Exp/Square.exp',300000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

%partitioning
md.qmu.numberofpartitions=md.mesh.numberofvertices;
md=partitioner(md,'package','linear','npart',md.qmu.numberofpartitions);
md.qmu.partition=md.qmu.partition-1;
md.qmu.isdakota=1;


%Dakota options

%dakota version
version=IssmConfig('_DAKOTA_VERSION_'); version=version(1:3); version=str2num(version);

%variables
md.qmu.variables.rho_ice=normal_uncertain('MaterialsRhoIce',md.materials.rho_ice,0.01);
md.qmu.variables.drag_coefficient=normal_uncertain('scaled_FrictionCoefficient',1,0.01);

%responses
md.qmu.responses.MaxVel=response_function('MaxVel',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);

%method
md.qmu.method     =dakota_method('nond_l');

%parameters
md.qmu.params.direct=true;
md.qmu.params.interval_type='forward';

if version>=6,
	md.qmu.params.analysis_driver='matlab';
	md.qmu.params.evaluation_scheduling='master';
	md.qmu.params.processors_per_evaluation=2;
else
	md.qmu.params.analysis_driver='stressbalance';
	md.qmu.params.evaluation_concurrency=1;
end


%imperative! 
md.stressbalance.reltol=10^-5; %tighten for qmu analysese

%solve
md=solve(md,'Stressbalance','overwrite','y');

%Fields and tolerances to track changes
md.qmu.results=md.results.dakota;
md.results.dakota.importancefactors=importancefactors(md,'scaled_FrictionCoefficient','MaxVel')';
field_names     ={'importancefactors'};
field_tolerances={1e-10};
field_values={...
         md.results.dakota.importancefactors,...
	};
