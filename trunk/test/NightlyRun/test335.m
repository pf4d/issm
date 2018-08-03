%Test Name: SquareSheetExtrudedHydrologyDCTwoLayers
md=triangle(model(),'../Exp/Square.exp',100000.);
md=setmask(md,'','');
md.transient=deactivateall(md.transient);
md.transient.ishydrology=1;
md.transient.issmb=1;
md=parameterize(md,'../Par/IceCube.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',1);
md.hydrology=(hydrologydc);
md.hydrology=initialize(md.hydrology,md);
md.hydrology.isefficientlayer=1;
md.hydrology.sedimentlimit_flag=1;
md.hydrology.sedimentlimit=800.0;
md.initialization.sediment_head=0.0*ones(md.mesh.numberofvertices,1);
md.hydrology.spcsediment_head=NaN*ones(md.mesh.numberofvertices,1);
md.basalforcings.groundedice_melting_rate = 2.0*ones(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate = 0.0*ones(md.mesh.numberofvertices,1);
md.hydrology.sediment_transmitivity=3*ones(md.mesh.numberofvertices,1);

md.initialization.epl_head=0.0*ones(md.mesh.numberofvertices,1);
md.initialization.epl_thickness=1.0*ones(md.mesh.numberofvertices,1);
md.hydrology.spcepl_head=NaN*ones(md.mesh.numberofvertices,1);
md.hydrology.mask_eplactive_node=0*ones(md.mesh.numberofvertices,1);
md.hydrology.epl_conductivity=30;
md.hydrology.epl_initial_thickness=1;
md.hydrology.epl_colapse_thickness=1.0e-3;
md.hydrology.epl_thick_comp=1;
md.hydrology.epl_max_thickness=1;
md.timestepping.time_step=0.2;
md.timestepping.final_time=2.0;

%md.verbose.solution=1;
md=extrude(md,3,1.1);
md=solve(md,'Transient');

%store=md.constants.g*md.hydrology.sediment_porosity* ...
%			md.materials.rho_freshwater*((md.hydrology.sediment_compressibility/md.hydrology.sediment_porosity)+md.hydrology.water_compressibility)

%sed=ones(1,size(md.results.HydrologySolution,2));
%epl=ones(1,size(md.results.HydrologySolution,2));
%res=ones(1,size(md.results.HydrologySolution,2));
%input=ones(1,size(md.results.HydrologySolution,2));
%for i= 1:size(md.results.HydrologySolution,2)
%	sed(i)=mean(md.results.HydrologySolution(i).SedimentHead);
%	res(i)=mean(md.results.HydrologySolution(i).SedimentHeadResidual);
%	epl(i)=mean(md.results.HydrologySolution(i).EplHead);
%	input(i)=2.0*(i*0.2);
%end

%Fields and tolerances to track changes
field_names     ={'SedimentWaterHead1','EplWaterHead1','SedimentHeadResidual1',...
						'SedimentWaterHead4','EplWaterHead4','SedimentHeadResidual4',...
						'SedimentWaterHead5','EplWaterHead5','SedimentHeadResidual5',...
						'SedimentWaterHead9','EplWaterHead9','SedimentHeadResidual9'};
field_tolerances={1e-13, 1e-13, 1e-13,...
						1e-13, 1e-13, 1e-13,...
						1e-13, 5e-12, 2e-11,...
						1e-13, 5e-12, 2e-11};
field_values={md.results.TransientSolution(1).SedimentHead, ...
							md.results.TransientSolution(1).EplHead,...
							md.results.TransientSolution(1).SedimentHeadResidual,...
							md.results.TransientSolution(4).SedimentHead,...
							md.results.TransientSolution(4).EplHead,...
							md.results.TransientSolution(4).SedimentHeadResidual, ...
							md.results.TransientSolution(5).SedimentHead,...
							md.results.TransientSolution(5).EplHead,...
							md.results.TransientSolution(5).SedimentHeadResidual, ...
							md.results.TransientSolution(9).SedimentHead,...
							md.results.TransientSolution(9).EplHead,...
							md.results.TransientSolution(9).SedimentHeadResidual};
