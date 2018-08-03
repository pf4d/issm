function exp2shp(expfilename,shpfilename,geometry)
%SHPWRITE - write a shape file from a contour structure
%
%   Usage:
%      exp2shp(expfilename,shpfilename,geometry)
%
%   Example:
%      exp2shp('domainoutline.exp','domainoutline.shp')
%      exp2shp('domainoutline.exp','domainoutline.shp','Polygon')
%      exp2shp('massfluxgate.exp','massfluxgate.shp','Line')
%
%   See also SHPREAD, SHPWRITE, SHP2EXP

%check file extensions
[pathstr,name,ext] = fileparts(shpfilename);
if ~strcmp(ext,'.shp'),
	error(['Shapefile ' shpfilename ' does not have an extension .shp']);
end

[pathstr,name,ext] = fileparts(expfilename);
if ~strcmp(ext,'.exp'),
	error(['Exp file ' expfilename ' does not have an extension .exp']);
end

shp=expread(expfilename);

%initialize number of profile
count=0;

contours=struct([]);
for i=1:length(shp),
	if nargin < 3
		if length(shp(1).x) == 1
			geometry = 'Point';
		elseif length(shp(1).x) < 3
			geometry = 'Line';
		else
			geometry = 'Polygon';
		end
	end
	contours(i).Geometry=geometry;
	contours(i).id=i;
	contours(i).X=shp(i).x;
	contours(i).Y=shp(i).y;
end
	
shapewrite(contours,shpfilename);
