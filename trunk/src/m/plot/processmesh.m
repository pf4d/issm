function [x y z elements is2d isplanet]=processmesh(md,data,options)
%PROCESSMESH - process mesh to be plotted
%
%   Usage:
%      [x y z elements is2d]=processmesh(md,data,options)
%
%   See also: PLOTMODEL, PROCESSDATA

%some checks
if md.mesh.numberofvertices==0,
	error('plot error message: mesh is empty')
end
if md.mesh.numberofvertices==md.mesh.numberofelements
	error(['plot error message: the number of elements is the same as the number of nodes...']);
end

%special case for mesg 2dvertical
if strcmp(domaintype(md.mesh),'2Dvertical'),
	[x y z elements is2d isplanet] = processmesh(md.mesh,options);
	return;
end

%special case for mesh 3dsurface
if strcmp(domaintype(md.mesh),'3Dsurface'),
	[x y z elements is2d isplanet] = processmesh(md.mesh,options);
	return;
end

if ~strcmpi(getfieldvalue(options,'coord','xy'),'latlon'),
	x=md.mesh.x;
	if isprop(md.mesh,'x2d'), x2d=md.mesh.x2d; end
	y=md.mesh.y;
	if isprop(md.mesh,'y2d'), y2d=md.mesh.y2d; end
else
	x=md.mesh.long;
	y=md.mesh.lat;
end

if isprop(md.mesh,'z'),
	z=md.mesh.z;
else
	z=zeros(size(x));
end
z=getfieldvalue(options,'z',z);
if ischar(z),
	z=md.(z);
end

if isprop(md.mesh,'elements2d'), elements2d=md.mesh.elements2d; end
elements=md.mesh.elements;

%is it a 2d plot?
if md.mesh.dimension()==2,
	is2d=1;
else
	if getfieldvalue(options,'layer',0)>=1,
		is2d=1;
	else
		is2d=0;
	end
end

%layer projection? 
if getfieldvalue(options,'layer',0)>=1,
	if strcmpi(getfieldvalue(options,'coord','xy'),'latlon'),
		error('processmesh error message: cannot work with 3D meshes for now');
	end
	%we modify the mesh temporarily to a 2d mesh from which the 3d mesh was extruded. 
	x=x2d;
	y=y2d;
	z=zeros(size(x2d));
	elements=elements2d;
end

%units
if exist(options,'unit'),
	unit=getfieldvalue(options,'unit');
	x=x*unit;
	y=y*unit;
	z=z*unit;
end

if isa(md,'planet'),
	isplanet=1;
else
	isplanet=0;
end
