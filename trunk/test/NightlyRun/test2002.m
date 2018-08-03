%Test Name: EarthSlr

%mesh earth: 
md=model; 
md.mesh=gmshplanet('radius',6.371012*10^3,'resolution',700.); %500 km resolution mesh

%parameterize slr solution:
%slr loading:  {{{
md.slr.deltathickness=zeros(md.mesh.numberofelements,1);
md.slr.sealevel=zeros(md.mesh.numberofvertices,1);
%antarctica
late=sum(md.mesh.lat(md.mesh.elements),2)/3;
longe=sum(md.mesh.long(md.mesh.elements),2)/3;
pos=find(late <-80);
md.slr.deltathickness(pos)=-100;
%greenland 
pos=find(late > 70 &  late < 80 & longe>-60 & longe<-30);
md.slr.deltathickness(pos)=-100;

%elastic loading from love numbers: 
nlov=101;
md.slr.love_h = love_numbers('h','CM'); md.slr.love_h(nlov+1:end)=[];
md.slr.love_k = love_numbers('k','CM'); md.slr.love_k(nlov+1:end)=[];
md.slr.love_l = love_numbers('l','CM'); md.slr.love_l(nlov+1:end)=[];

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

md.slr.ocean_area_scaling=0; 

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
md.miscellaneous.name='test2002';

%Solution parameters
md.slr.reltol=NaN;
md.slr.abstol=1e-3;

%eustatic run: 
md.slr.rigid=0; md.slr.elastic=0;
md=solve(md,'Sealevelrise');
Seustatic=md.results.SealevelriseSolution.Sealevel;

%eustatic + rigid run: 
md.slr.rigid=1; md.slr.elastic=0;
md=solve(md,'Sealevelrise');
Srigid=md.results.SealevelriseSolution.Sealevel;

%eustatic + rigid + elastic run: 
md.slr.rigid=1; md.slr.elastic=1;
md=solve(md,'Sealevelrise');
Selastic=md.results.SealevelriseSolution.Sealevel;

%eustatic + rigid + elastic + rotation run: 
md.slr.rigid=1; md.slr.elastic=1; md.slr.rotation=1;
md=solve(md,'Sealevelrise');
Srotation=md.results.SealevelriseSolution.Sealevel;

%Fields and tolerances to track changes
field_names     ={'Eustatic','Rigid','Elastic','Rotation'};
field_tolerances={1e-13,1e-13,1e-13,1e-13};
field_values={Seustatic,Srigid,Selastic,Srotation};
