function areas=GetAreas3DTria(index,x,y,z,varargin)
%GETAREAS3DTRIA - compute areas of triangles with 3D coordinates 
%
%   compute areas of trianguls with 3D coordinates 
%
%   Usage:
%      areas  =GetAreas3DTria(index,x,y,z);
%
%   Examples:
%      areas  =GetAreas3DTria(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z);

%get number of elements and number of nodes
nels=size(index,1); 
nods=length(x);  

%some checks
if nargout~=1 | (nargin~=3 & nargin~=4),
	help GetAreas3DTria
	error('GetAreas error message: bad usage')
end
if ((length(y)~=nods) | (nargin==4 & length(z)~=nods)),
	error('GetAreas3DTria error message: x,y and z do not have the same length')
end
if max(index(:))>nods,
	error(['GetAreas3DTria error message: index should not have values above ' num2str(nods) ])
end
if (nargin==4 & size(index,2)~=3),
	error('GetAreas3DTria error message: index should have 3 columns for 2d meshes.')
end

%initialization
areas=zeros(nels,1);
x1=x(index(:,1)); x2=x(index(:,2)); x3=x(index(:,3));
y1=y(index(:,1)); y2=y(index(:,2)); y3=y(index(:,3));
z1=z(index(:,1)); z2=z(index(:,2)); z3=z(index(:,3));

%compute the volume of each element
if nargin==4,
   % area of triangles with 3D coordinats
   for j=1:nels
		m1=[x1(j) x2(j) x3(j); y1(j) y2(j) y3(j); 1 1 1];
	   m2=[y1(j) y2(j) y3(j); z1(j) z2(j) z3(j); 1 1 1];
	   m3=[z1(j) z2(j) z3(j); x1(j) x2(j) x3(j); 1 1 1];
      areas(j)=sqrt(det(m1)^2 + det(m2)^2 + det(m3)^2)/2;
   end
end 


