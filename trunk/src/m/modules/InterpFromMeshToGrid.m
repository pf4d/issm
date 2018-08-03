function [x_m,y_m,griddata] = InterpFromMeshToGrid(index,x,y,data,xmin,ymax,xposting,yposting,nlines,ncols,default_value);
%INTERPFROMMESHTOGRID - Interpolation of a data defined on a mesh onto a grid
%
%   This function is a multi-threaded mex file that interpolates a field defined on a triangular
%   mesh onto a regular grid.
%
%   index,x,y:	delaunay triangulation defining the mesh
%   meshdata:	vertex values of data to be interpolated
%
%   xmin,ymax,posting,nlines,ncols: parameters that define the grid
%   default_value:	value of points located out of the mesh

% Check usage
if nargin~=11
	help InterpFromMeshToGrid
	error('Wrong usage (see above)');
end

% Call mex module
[x_m,y_m,griddata] = InterpFromMeshToGrid_matlab(index,x,y,data,xmin,ymax,xposting,yposting,nlines,ncols,default_value);
