%Test Name: Esa2Dsurface
%Elastostatic adjustment for Wahr problem (Wahr et al., 2013, Fig. 1) 

%mesh earth: 
md=model;
md=roundmesh(md,50000,2000);     % radius and element size (meters)

% define load 
md.esa.deltathickness=zeros(md.mesh.numberofelements,1);
disc_radius=20; % km
index=md.mesh.elements;
x_element=mean(md.mesh.x(index),2);
y_element=mean(md.mesh.y(index),2);
rad_dist=sqrt(x_element.^2+y_element.^2)/1000;  % radial distance in km
md.esa.deltathickness(rad_dist<=disc_radius)=-1;   % 1 m water withdrawl

%love numbers:
nlov=10001;	% horizontal displacements do not work for low degree truncation, e.g., 101
md.esa.love_h = love_numbers('h'); md.esa.love_h(nlov+1:end)=[];
md.esa.love_l = love_numbers('l'); md.esa.love_l(nlov+1:end)=[];

%mask:  {{{
	%make sure wherever there is an ice load, that the mask is set to ice: 
	md.mask.ice_levelset=ones(md.mesh.numberofvertices,1);
	pos=find(md.esa.deltathickness); md.mask.ice_levelset(md.mesh.elements(pos,:))=-1;

	%is ice grounded? 
	md.mask.groundedice_levelset=-ones(md.mesh.numberofvertices,1);
	pos=find(md.mask.ice_levelset<=0); md.mask.groundedice_levelset(pos)=1; 

% }}}
%geometry:  {{{
	di=md.materials.rho_ice/md.materials.rho_water;
	md.geometry.thickness=ones(md.mesh.numberofvertices,1);
	md.geometry.surface=(1-di)*zeros(md.mesh.numberofvertices,1);
	md.geometry.base=md.geometry.surface-md.geometry.thickness;
	md.geometry.bed=md.geometry.base;
% }}}
%materials:  {{{
	md.initialization.temperature=273.25*ones(md.mesh.numberofvertices,1);
	md.materials.rheology_B=paterson(md.initialization.temperature);
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
% }}}
%Miscellaneous: {{{
	md.miscellaneous.name='test2110';
% }}}

md.esa.degacc=0.01;

% solve esa 
md.esa.requested_outputs = {'EsaUmotion','EsaNmotion','EsaEmotion',...
	'EsaStrainratexx','EsaStrainratexy','EsaStrainrateyy','EsaRotationrate'};
md.cluster=generic('name',oshostname(),'np',3); 
md.verbose=verbose('111111111');
md=solve(md,'Esa');

%Fields and tolerances to track changes
field_names     ={'EsaUmotion','EsaNmotion','EsaEmotion',...
	'EsaStrainratexx','EsaStrainratexy','EsaStrainrateyy','EsaRotationrate'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.EsaSolution.EsaUmotion),...
	(md.results.EsaSolution.EsaNmotion),...
	(md.results.EsaSolution.EsaEmotion),...
	(md.results.EsaSolution.EsaStrainratexx),...
	(md.results.EsaSolution.EsaStrainratexy),...
	(md.results.EsaSolution.EsaStrainrateyy),...
	(md.results.EsaSolution.EsaRotationrate),...
	};

