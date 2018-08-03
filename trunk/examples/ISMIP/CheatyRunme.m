%wich steps to perform steps are from 1 to 8 
%step 7 is specific to ISMIPA
%step 8 is specific to ISMIPF

steps=[1];

% parameter file to be used, choose between CheatyIsmipA.par or CheatyIsmipF.par
ParamFile='CheatyIsmipF.par'

%Run Steps

%Mesh Generation #1
if any(steps==1)

	%initialize md as a new model #help model
	%->
	md=model();
	% generate a squaremesh #help squaremesh
	% Side is 80 km long with 20 points
	%->
	if(ParamFile=='CheatyIsmipA.par'),
		md=squaremesh(md,80000,80000,20,20);
	elseif(ParamFile=='CheatyIsmipF.par'),
		md=squaremesh(md,100000,100000,30,30);
  end
	% plot the given mesh #plotdoc
	%->
	plotmodel(md,'data','mesh')
	% save the given model
	%->
	save ./Models/ISMIP.Mesh_generation md;
end

%Masks #2
if any(steps==2)

	% load the preceding step #help loadmodel
	% path is given by the organizer with the name of the given step
	%->
	md = loadmodel('./Models/ISMIP.Mesh_generation');
	% set the mask #help setmask
	% all MISMIP nodes are grounded
	%->
	md=setmask(md,'','');
	% plot the given mask #md.mask to locate the field
	%->
	plotmodel(md,'data',md.mask.groundedice_levelset);
	% save the given model
	%->
	save ./Models/ISMIP.SetMask md;
end

%Parameterization #3
if any(steps==3)

	% load the preceding step #help loadmodel
	% path is given by the organizer with the name of the given step
	%->
	md = loadmodel('./Models/ISMIP.SetMask');
	% parametrize the model # help parameterize
	% you will need to fil-up the parameter file defined by the
	% ParamFile variable
	%->
	md=parameterize(md,ParamFile);
	% save the given model
	%->
	save ./Models/ISMIP.Parameterization md;
end

%Extrusion #4
if any(steps==4)
	
	% load the preceding step #help loadmodel
	% path is given by the organizer with the name of the given step
	%->
	md = loadmodel('./Models/ISMIP.Parameterization');
	% vertically extrude the preceding mesh #help extrude
	% only 5 layers exponent 1
	%->
	md=extrude(md,9,1);
	% plot the 3D geometry #plotdoc
	%->
	plotmodel(md,'data',md.geometry.base)
	% save the given model
	%->
	save ./Models/ISMIP.Extrusion md;
end

%Set the flow computing method #5
if any(steps==5)

	% load the preceding step #help loadmodel
	% path is given by the organizer with the name of the given step
	%->
	md = loadmodel('./Models/ISMIP.Extrusion');
	% set the approximation for the flow computation #help setflowequation
	% We will be using the Higher Order Model (HO)
	%->
	md=setflowequation(md,'HO','all');
	% save the given model
	%->
	save ./Models/ISMIP.SetFlow md;
end

%Set Boundary Conditions #6
if any(steps==6)

	% load the preceding step #help loadmodel
	% path is given by the organizer with the name of the given step
	%->
	md = loadmodel('./Models/ISMIP.SetFlow');
	% dirichlet boundary condition are known as SPCs
	% ice frozen to the base, no velocity	#md.stressbalance
	% SPCs are initialized at NaN one value per vertex
	%->
	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	%->
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	%->
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
	% extract the nodenumbers at the base #md.mesh.vertexonbase
	%->
	basalnodes=find(md.mesh.vertexonbase);
	% set the sliding to zero on the bed
	%->
	md.stressbalance.spcvx(basalnodes)=0.0;
	%->
	md.stressbalance.spcvy(basalnodes)=0.0;
	% periodic boundaries have to be fixed on the sides
	% create tabs with the side of the domain
	% for x
	% create maxX #help find
	%->
	maxX=find(md.mesh.x==max(md.mesh.x));
	% create minX
	%->
	minX=find(md.mesh.x==min(md.mesh.x));
	% for y, max X and minX should be excluded
	% create maxY
	%->
	maxY=find(md.mesh.y==max(md.mesh.y) & md.mesh.x~=max(md.mesh.x) & md.mesh.x~=min(md.mesh.x));
	% create minY
	%->
	minY=find(md.mesh.y==min(md.mesh.y) & md.mesh.x~=max(md.mesh.x) & md.mesh.x~=min(md.mesh.x));
	% set the node that should be paired together
	% #md.stressbalance.vertex_pairing
	%->
	md.stressbalance.vertex_pairing=[minX,maxX;minY,maxY];
	if (ParamFile=='CheatyIsmipF.par')
		% if we are dealing with IsmipF the solution is in
		% masstransport
		md.masstransport.vertex_pairing=md.stressbalance.vertex_pairing;
  end
	% save the given model
	%->
	save ./Models/ISMIP.BoundaryCondition md;
end

%Solving #7
if any(steps==7)
	% load the preceding step #help loadmodel
	% path is given by the organizer with the name of the given step
	%->
	md = loadmodel('./Models/ISMIP.BoundaryCondition');
	% Set cluster #md.cluster
	% generic parameters #help generic
	% set only the name and number of process
	%->
	md.cluster=generic('name',oshostname(),'np',2);
	% Set which control message you want to see #help verbose
	%->
	md.verbose=verbose('convergence',true);
	% Solve #help solve
	% we are solving a StressBalanc
	%->
	md=solve(md,'Stressbalance');
	% save the given model
	%->
	save ./Models/ISMIP.StressBalance md;
	% plot the surface velocities #plotdoc
	%->
	plotmodel(md,'data',md.results.StressbalanceSolution.Vel)
end

%Solving #8
if any(steps==8)
	% load the preceding step #help loadmodel
	% path is given by the organizer with the name of the given step
	%->
	md = loadmodel('./Models/ISMIP.BoundaryCondition');
	% Set cluster #md.cluster
	% generic parameters #help generic
	% set only the name and number of process
	%->
	md.cluster=generic('name',oshostname(),'np',2);
	% Set which control message you want to see #help verbose
	%->
	md.verbose=verbose('convergence',true);
	% set the transient model to ignore the thermal model
	% #md.transient 
	%->
	md.transient.isthermal=0;
	% define the timestepping scheme
	% everything here should be provided in years #md.timestepping
	% give the length of the time_step (4 years)
	%->
	md.timestepping.time_step=4;
	% give final_time (20*4 years time_steps)
	%->
	md.timestepping.final_time=4*20;
	% Solve #help solve
	% we are solving a TransientSolution
	%->
	md=solve(md,'Transient');
	% save the given model
	%->
	save ./Models/ISMIP.Transient md;
	% plot the surface velocities #plotdoc
	%->
	plotmodel(md,'data',md.results.TransientSolution(20).Vel)
end
