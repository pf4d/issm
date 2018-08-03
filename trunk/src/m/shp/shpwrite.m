function shpwrite(shp,filename)
%SHPWRITE - write a shape file from a contour structure
%
%   Usage:
%      shpwrite(shp,filename)
%
%   Example:
%      shpwrite(shp,'domainoutline.shp')
%
%   See also SHPREAD


%initialize number of profile
count=0;

contours=struct([]);
for i=1:length(shp),
	if strcmpi(shp(i).Geometry,'Point'),
		contours(i).Geometry='Point';
	else strcmpi(shp(i).Geometry,'Polygon'),
		contours(i).Geometry='Polygon';
	end
	contours(i).id=i;
	contours(i).X=shp(i).x;
	contours(i).Y=shp(i).y;
end
	
shapewrite(contours,filename);
end
