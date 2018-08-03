%Test Name: SquareShelfTranIspddIsdeltaSSA2d
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');

%md.verbose=verbose('all');

% Use of ispdd and isdelta18o methods
md.smb = SMBpdd();
md.smb.isdelta18o=1;
md.smb.ismungsm=0;

%md.smb.precipitation(1:md.mesh.numberofvertices,1:12)=0;
%md.smb.monthlytemperatures(1:md.mesh.numberofvertices,1:12)=273;

% Add temperature, precipitation and delta18o needed to measure the surface mass balance
%  creating delta18o
load '../Data/delta18o.data'
md.smb.delta18o=delta18o;
% creating delta18oSurface
md.smb.delta18o_surface(1,1:(length(delta18o))) = 0;
md.smb.delta18o_surface(2,:) = delta18o(2,:);

% creating Present day and lgm temperatures
% Same temperature over the all region:
tmonth(1:12)=238.15+20.;
for imonth=0:11
    md.smb.temperatures_presentday(1:md.mesh.numberofvertices,imonth+1)=tmonth(imonth+1);
    md.smb.temperatures_lgm(1:md.mesh.numberofvertices,imonth+1)=tmonth(imonth+1)-20.;
    % Time for the last line:
    md.smb.temperatures_presentday(md.mesh.numberofvertices+1,imonth+1)=((imonth+1)/12);
    md.smb.temperatures_lgm(md.mesh.numberofvertices+1,imonth+1)=((imonth+1)/12);
end

% creating initialization and spc temperatures initialization and
% spc
md.thermal.spctemperature=mean(md.smb.temperatures_lgm(1:md.mesh.numberofvertices,1:12),2); %-10*ones(md.mesh.numberofvertices,1);
md.thermal.spctemperature=repmat(md.thermal.spctemperature,1,md.timestepping.final_time/md.timestepping.time_step);
itemp=0:md.timestepping.time_step:md.timestepping.final_time-md.timestepping.time_step;
md.thermal.spctemperature(md.mesh.numberofvertices+1,:)=itemp;

md.initialization.temperature=md.smb.temperatures_lgm(1:md.mesh.numberofvertices,1); %*ones(md.mesh.numberofvertices,1);
md.smb = initialize(md.smb,md);

% creating precipitation
for imonth=0:11
    md.smb.precipitations_presentday(1:md.mesh.numberofvertices,imonth+1)=-0.4*10^(-6)*md.mesh.y+0.5;
    md.smb.precipitations_lgm(1:md.mesh.numberofvertices,imonth+1)=-0.4*10^(-6)*md.mesh.y+0.5;
    % Time for the last line:
    md.smb.precipitations_presentday(md.mesh.numberofvertices+1,imonth+1)=((imonth+1)/12);
    md.smb.precipitations_lgm(md.mesh.numberofvertices+1,imonth+1)=((imonth+1)/12);
end

% Interpolation factors
md.smb.Tdiff(1,1:md.timestepping.final_time)=0.5;
md.smb.sealev(1,1:md.timestepping.final_time)=0.5;
% Year of each data point
md.smb.Tdiff(2,1:md.timestepping.final_time)=1:1:md.timestepping.final_time;
md.smb.sealev(2,1:md.timestepping.final_time)=1:1:md.timestepping.final_time;

% time steps and resolution
md.timestepping.time_step=20;
md.settings.output_frequency=1;
md.timestepping.final_time=60;

% 
md.transient.requested_outputs={'default','SmbMonthlytemperatures'};
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',1); % 3 for the cluster
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'Vx1','Vy1','Vel1','Pressure1','Bed1','Surface1','Thickness1','SmbMonthlytemperatures1','SmbMassBalance1',...
	   'Vx2','Vy2','Vel2','Pressure2','Bed2','Surface2','Thickness2','SmbMonthlytemperatures2','SmbMassBalance2',...
	   'Vx3','Vy3','Vel3','Pressure3','Bed3','Surface3','Thickness3','SmbMonthlytemperatures3','SmbMassBalance3'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,...
	1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,...
	1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.TransientSolution(1).Vx),...
	(md.results.TransientSolution(1).Vy),...
	(md.results.TransientSolution(1).Vel),...
	(md.results.TransientSolution(1).Pressure),...
	(md.results.TransientSolution(1).Base),...
	(md.results.TransientSolution(1).Surface),...
	(md.results.TransientSolution(1).Thickness),...
	(md.results.TransientSolution(1).SmbMonthlytemperatures),...
	(md.results.TransientSolution(1).SmbMassBalance),...
	(md.results.TransientSolution(2).Vx),...
	(md.results.TransientSolution(2).Vy),...
	(md.results.TransientSolution(2).Vel),...
	(md.results.TransientSolution(2).Pressure),...
	(md.results.TransientSolution(2).Base),...
	(md.results.TransientSolution(2).Surface),...
	(md.results.TransientSolution(2).Thickness),...
	(md.results.TransientSolution(2).SmbMonthlytemperatures),...
	(md.results.TransientSolution(2).SmbMassBalance),...
	(md.results.TransientSolution(3).Vx),...
	(md.results.TransientSolution(3).Vy),...
	(md.results.TransientSolution(3).Vel),...
	(md.results.TransientSolution(3).Pressure),...
	(md.results.TransientSolution(3).Base),...
	(md.results.TransientSolution(3).Surface),...
	(md.results.TransientSolution(3).Thickness),...
	(md.results.TransientSolution(3).SmbMonthlytemperatures),...
	(md.results.TransientSolution(3).SmbMassBalance),...
	};
