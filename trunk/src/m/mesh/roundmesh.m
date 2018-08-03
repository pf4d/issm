function md=roundmesh(md,radius,resolution)
%ROUNDMESH - create an unstructured round mesh 
%
%   This script will generate a structured round mesh
%   - radius     : specifies the radius of the circle in meters
%   - resolution : specifies the resolution in meters
%
%   Usage:
%      md=roundmesh(md,radius,resolution)

%First we have to create the domain outline 

%Get number of points on the circle
pointsonedge=floor((2.*pi*radius) / resolution);

%Calculate the cartesians coordinates of the points
x_list=ones(pointsonedge,1); y_list=ones(pointsonedge,1);
theta=(0.:2.*pi/pointsonedge:2.*pi*(1.-1./pointsonedge))';
x_list=roundsigfig(radius*x_list.*cos(theta),12);
y_list=roundsigfig(radius*y_list.*sin(theta),12);
A=struct('x',x_list,'y',y_list,'density',1.);
expwrite(A,'RoundDomainOutline.exp');

%Call Bamg
md=triangle(md,'RoundDomainOutline.exp',resolution);
%md=bamg(md,'domain','RoundDomainOutline.exp','hmin',resolution);

%move the closest node to the center
[mini pos]=min(md.mesh.x.^2+md.mesh.y.^2);
md.mesh.x(pos)=0.;
md.mesh.y(pos)=0.;

%delete domain
delete('RoundDomainOutline.exp')
end

function x=roundsigfig(x,n)

digits=ceil(log10(abs(x)));
x=x./10.^digits;
x=round(x.*10.^n)./10.^n;
x=x.*10.^digits;

pos=find(isnan(x));
x(pos)=0.;

end
