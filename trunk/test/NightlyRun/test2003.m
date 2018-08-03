%Test Name: EarthSlr_rotationalFeedback

%mesh earth: 
md=model; 
md.mesh=gmshplanet('radius',6.371012*10^3,'resolution',1000.); %500 km resolution mesh

%parameterize slr solution:
%slr loading:  {{{
md.slr.deltathickness=zeros(md.mesh.numberofelements,1);
md.slr.sealevel=zeros(md.mesh.numberofvertices,1);
%antarctica
late=sum(md.mesh.lat(md.mesh.elements),2)/3;
longe=sum(md.mesh.long(md.mesh.elements),2)/3;
pos=find(late <-75 & longe <0);
md.slr.deltathickness(pos)=-1;

%elastic loading from love numbers: 
nlov=1001;
md.slr.love_h = love_numbers('h'); md.slr.love_h(nlov+1:end)=[];
md.slr.love_k = love_numbers('k'); md.slr.love_k(nlov+1:end)=[];
md.slr.love_l = love_numbers('l'); md.slr.love_l(nlov+1:end)=[];

%}}}
%mask:  {{{
md.mask=maskpsl(); % use maskpsl class (instead of mask) to store the ocean function as a ocean_levelset 
mask=gmtmask(md.mesh.lat,md.mesh.long); 

icemask=ones(md.mesh.numberofvertices,1);
pos=find(mask==0);  icemask(pos)=-1;
pos=find(sum(mask(md.mesh.elements),2)<3);   icemask(md.mesh.elements(pos,:))=-1;
md.mask.ice_levelset=icemask;
md.mask.ocean_levelset=zeros(md.mesh.numberofvertices,1);
pos=find(md.mask.ice_levelset==1); md.mask.ocean_levelset(pos)=1;

%make sure that the ice level set is all inclusive:
md.mask.land_levelset=zeros(md.mesh.numberofvertices,1);
md.mask.groundedice_levelset=-ones(md.mesh.numberofvertices,1); 

%make sure wherever there is an ice load, that the mask is set to ice: 
pos=find(md.slr.deltathickness); md.mask.ice_levelset(md.mesh.elements(pos,:))=-1;
% }}}

% use model representation of ocea area (not the ture area) 
md.slr.ocean_area_scaling = 0; 

%geometry
di=md.materials.rho_ice/md.materials.rho_water;
md.geometry.thickness=ones(md.mesh.numberofvertices,1);
md.geometry.surface=(1-di)*zeros(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.surface-md.geometry.thickness;
md.geometry.bed=md.geometry.base;

%materials
md.initialization.temperature=273.25*ones(md.mesh.numberofvertices,1);
md.materials.rheology_B=paterson(md.initialization.temperature);
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

%Miscellaneous
md.miscellaneous.name='test2003';

%Solution parameters
md.slr.reltol=NaN;
md.slr.abstol=1e-3;

%eustatic + rigid + elastic run: 
md.slr.rigid=1; md.slr.elastic=1; md.slr.rotation=0; 
md.cluster=generic('name',oshostname(),'np',3);
%md.verbose=verbose('111111111');
md=solve(md,'Sealevelrise');
SnoRotation=md.results.SealevelriseSolution.Sealevel;

%eustatic + rigid + elastic + rotation run: 
md.slr.rigid=1; md.slr.elastic=1; md.slr.rotation=1;
md.cluster=generic('name',oshostname(),'np',3);
%md.verbose=verbose('111111111');
md=solve(md,'Sealevelrise');
SRotation=md.results.SealevelriseSolution.Sealevel;

%Fields and tolerances to track changes
field_names     ={'noRotation','Rotation'};
field_tolerances={1e-13,1e-13};
field_values={SnoRotation,SRotation};
