%Test Name: SquareShelfSMBGemb
md=triangle(model(),'../Exp/Square.exp',200000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=setflowequation(md,'SSA','all');
md.materials.rho_ice=910;
md.cluster=generic('name',oshostname(),'np',3);

% Use of Gemb method for SMB computation
md.smb = SMBgemb(md.mesh,md.geometry);

%load hourly surface forcing date from 1979 to 2009:
inputs=load('../Data/gemb_input.mat');

%setup the inputs: 
md.smb.Ta=[repmat(inputs.Ta0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.V=[repmat(inputs.V0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.dswrf=[repmat(inputs.dsw0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.dlwrf=[repmat(inputs.dlw0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.P=[repmat(inputs.P0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.eAir=[repmat(inputs.eAir0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.pAir=[repmat(inputs.pAir0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.Vz=repmat(inputs.LP.Vz,md.mesh.numberofelements,1);
md.smb.Tz=repmat(inputs.LP.Tz,md.mesh.numberofelements,1);
md.smb.Tmean=repmat(inputs.LP.Tmean,md.mesh.numberofelements,1);
md.smb.C=repmat(inputs.LP.C,md.mesh.numberofelements,1);

%smb settings
md.smb.requested_outputs={'SmbDz','SmbT','SmbD','SmbRe','SmbGdn','SmbGsp','SmbEC','SmbA','SmbMassBalance'};

%only run smb core: 
md.transient.isstressbalance=0;
md.transient.ismasstransport=0;
md.transient.isthermal=0;

%time stepping: 
md.timestepping.start_time=1965;
md.timestepping.final_time=1966;
md.timestepping.time_step=1/365.0;
md.timestepping.interp_forcings=0;

%Run transient
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names      ={'SmbDz','SmbT' ,'SmbD' ,'SmbRe','SmbGdn','SmbGsp','SmbA' ,'SmbEC','SmbMassBalance'};
field_tolerances ={5e-4,5e-5,0.0006,0.0002,1e-5,0.0003,2e-5,2e-7,1e-7};

field_values={...
	(md.results.TransientSolution(end).SmbDz(1,1:240)),...
	(md.results.TransientSolution(end).SmbT(1,1:240)),...
	(md.results.TransientSolution(end).SmbD(1,1:240)),...
	(md.results.TransientSolution(end).SmbRe(1,1:240)),...
	(md.results.TransientSolution(end).SmbGdn(1,1:240)),...
	(md.results.TransientSolution(end).SmbGsp(1,1:240)),...
	(md.results.TransientSolution(end).SmbA(1,1:240)),...
	(md.results.TransientSolution(end).SmbEC(1)),...
	(md.results.TransientSolution(end).SmbMassBalance(1)),...
	};
