function [index,x,y,segments,segmentmarkers] = TriMesh(domainoutlinefilename,rifts,mesh_area);
%TRIMESH - Mesh a domain using an .exp file
%
%   Usage: 
%     [index,x,y,segments,segmentmarkers]=TriMesh(domainoutlinefilename,rifts,mesh_area);
%	      
%   index,x,y:	Defines a triangulation 
%   segments:	Array made of exterior segments to the mesh domain outline 
%   segmentmarkers:	Array flagging each segment
%
%   domainoutlinefilename:	Argus domain outline file
%   mesh_area:	Maximum area desired for any element of the resulting mesh

% Check usage
if nargin~=3 && nargout~=5
	help TriMesh
	error('Wrong usage (see above)');
end

% Call mex module
[index,x,y,segments,segmentmarkers]=TriMesh_matlab(domainoutlinefilename,rifts,mesh_area);
