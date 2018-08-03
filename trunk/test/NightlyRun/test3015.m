%Test Name: SquareShelfConstrainedMasstransp2dAdolcForwardDifference
%This test runs test3005 with autodiff on, and checks that 
%the value of the scalar forward difference match a step-wise differential

%First configure
md=triangle(model(),'../Exp/Square.exp',50000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',1);
md.masstransport.requested_outputs={'IceVolume'};
md.verbose=verbose('autodiff',true);
md.toolkits.DefaultAnalysis=issmgslsolver();

%setup autodiff parameters
index=1; %this is the scalar component we are checking against
md.autodiff.independents={...
	independent('name','md.geometry.thickness','type','vertex','nods',md.mesh.numberofvertices,'fos_forward_index',index)
	};

md.autodiff.dependents={...
	dependent('name','IceVolume','type','scalar')...
	};
md.autodiff.driver='fos_forward';

%parameters for the step-wise derivative
delta=0.001;
h1=md.geometry.thickness(index);
h0=h1*(1.-delta);
h2=h1*(1.+delta);
deltaH=(h2-h0);

%save model:
md2=md;

%evaluate derivative by forward and backward stepping 
%forward
md=md2;
md.autodiff.isautodiff=false;
md.geometry.thickness(index)=h0;
md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
md.geometry.surface=md.geometry.base+md.geometry.thickness;
md=SetIceShelfBC(md);

md=solve(md,'Masstransport');
V0=md.results.MasstransportSolution.IceVolume;

%backward
md=md2;
md.autodiff.isautodiff=false;
md.geometry.thickness(index)=h2;
md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
md.geometry.surface=md.geometry.base+md.geometry.thickness;
md=SetIceShelfBC(md);

md=solve(md,'Masstransport');
V2=md.results.MasstransportSolution.IceVolume;

%compute resulting derivative
dVdh_an=(V2-V0)/deltaH;

%evaluate derivative using ADOLC 
md=md2;
md.autodiff.isautodiff=true;
md.geometry.thickness(index)=h1;
md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
md.geometry.surface=md.geometry.base+md.geometry.thickness;
md=SetIceShelfBC(md);

md=solve(md,'Masstransport');
%retrieve directly
dVdh_ad=md.results.MasstransportSolution.AutodiffJacobian;

disp(sprintf('dV/dh: analytical:  %16.16g\n       using adolc:  %16.16g\n',dVdh_an,dVdh_ad));

%Fields and tolerances to track changes
field_names     ={'dV/dh-dV/dh0'};
field_tolerances={1e-13};
field_values={dVdh_ad-dVdh_an};
