function ShpRead(filename);
%	SHPREAD - Read shapefile
%	
%	   This module reads shapefiles and converts them to matlab/python structures
%	
%	   Usage:
%	      ShpRead(filename);
%	      filexp:      file name of exp file to be written
%	
%	   Examples:
%	      ShpRead('file.shp');

% Check usage
if nargin~=1
	help ShpRead
	error('Wrong usage: No file specified');
end

% Call mex module
ShpRead_matlab(filename);
