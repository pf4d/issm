%Test Name: EarthEsa
%Elastostatic adjustment for an elemental ice unloading 

%mesh earth: 
md=model; 
md.mesh=gmshplanet('radius',6.371012*10^3,'resolution',1000);

% define load 
md.esa.deltathickness=zeros(md.mesh.numberofelements,1);
pos=450;    md.esa.deltathickness(pos)=-100;   % this is the only "icy" element

%love numbers:
nlov=10001;
md.esa.love_h = love_numbers('h','CM'); md.esa.love_h(nlov+1:end)=[];
md.esa.love_l = love_numbers('l','CM'); md.esa.love_l(nlov+1:end)=[];

%mask:  {{{
	md.mask=maskpsl(); % use maskpsl class (instead of mask) to store the ocean function as a ocean_levelset 
	md.mask.ocean_levelset=gmtmask(md.mesh.lat,md.mesh.long); 

	%make sure wherever there is an ice load, that the mask is set to ice: 
	md.mask.ice_levelset=ones(md.mesh.numberofvertices,1);
	pos=find(md.esa.deltathickness); md.mask.ice_levelset(md.mesh.elements(pos,:))=-1;

	%is ice grounded? 
	md.mask.groundedice_levelset=-ones(md.mesh.numberofvertices,1);
	pos=find(md.mask.ice_levelset<=0); md.mask.groundedice_levelset(pos)=1; 

	%make sure ice domain is on the continent: 
	%pos=find(md.mask.ice_levelset<=0); md.mask.ocean_levelset(pos)=0; 

	%land mask 
	md.mask.land_levelset=1-md.mask.ocean_levelset; 

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
	md.miscellaneous.name='test2101';
% }}}

% solve esa 
md.esa.requested_outputs = {'EsaUmotion','EsaNmotion','EsaEmotion'};
md.cluster=generic('name',oshostname(),'np',3); 
md.verbose=verbose('111111111');
md=solve(md,'Esa');

%Fields and tolerances to track changes
field_names     ={'EsaUmotion','EsaNmotion','EsaEmotion'};
field_tolerances={1e-13,1e-13,1e-13};
field_values={...
	(md.results.EsaSolution.EsaUmotion),...
	(md.results.EsaSolution.EsaNmotion),...
	(md.results.EsaSolution.EsaEmotion),...
	};

