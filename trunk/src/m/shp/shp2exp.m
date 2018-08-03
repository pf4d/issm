function shp2exp(shpfilename,expfilename)
%SHP2EXP- transform shape file to Argus .exp file
%
%   Usage:
%      shp2exp(shpfilename,expfilename);
%
%   Example:
%      shp2exp('Domain.shp','Domain.exp');
%
%   See also EXPMASTER, EXPDOC

%check file extensions
[pathstr,name,ext] = fileparts(shpfilename);
if ~strcmp(ext,'.shp'),
	error(['Shapefile ' shpfilename ' does not have an extension .shp']);
end

[pathstr,name,ext] = fileparts(expfilename);
if ~strcmp(ext,'.exp'),
	error(['Exp file ' expfilename ' does not have an extension .exp']);
end

if ~exist(shpfilename,'file'),
	error(['Shapefile ' shpfilename ' does not exist']);
end
shp=shaperead(shpfilename);

expstruct=struct([]);
for i=1:length(shp),
	if strcmpi(shp(i).Geometry,'Polygon'),
		x=shp(i).X; y=shp(i).Y;
		ids=find(isnan(x));
		x(ids)=[]; y(ids)=[];
		expstruct(end+1).x=x;
		expstruct(end).y=y;
		expstruct(end).nods=length(x);
		expstruct(end).density=1;
		expstruct(end).closed=1;
		expstruct(end).name=num2str(shp(i).id);
	elseif strcmpi(shp(i).Geometry,'Point'),
		x=shp(i).X; y=shp(i).Y;
		expstruct(end+1).x=x;
		expstruct(end).y=y;
		expstruct(end).nods=length(x);
		expstruct(end).density=1;
		expstruct(end).closed=0;
		%exp(end).name=num2str(shp(i).id);
	elseif strcmpi(shp(i).Geometry,'Line'),
		x=shp(i).X; y=shp(i).Y;
		expstruct(end+1).x=x;
		expstruct(end).y=y;
		expstruct(end).nods=length(x);
		expstruct(end).density=1;
		expstruct(end).closed=0;
	end
end

expwrite(expstruct,expfilename);
