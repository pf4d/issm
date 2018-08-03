function md=squaremesh(md,Lx,Ly,nx,ny)
%SQUAREMESH - create a structured square mesh 
%
%   This script will generate a structured square mesh
%   Lx and Ly are the dimension of the domain (in meters)
%   nx anx ny are the number of nodes in the x and y direction
%   The coordinates x and y returned are in meters.
%
%   Usage:
%      [md]=squaremesh(md,Lx,Ly,nx,ny)

%get number of elements and number of nodes
nel=(nx-1)*(ny-1)*2;
nods=nx*ny;

%initialization
index=zeros(nel,3);
x=zeros(nx*ny,1);
y=zeros(nx*ny,1);

%create coordinates
for n=1:nx,
	for m=1:ny,
		x((n-1)*ny+m)=(n-1.);
		y((n-1)*ny+m)=(m-1.);
	end
end

%create index
for n=1:(nx-1)
	for m=1:(ny-1),
		A=(n-1)*ny+m;
		B=A+1;
		C=n*ny+m;
		D=C+1;
		index((n-1)*(ny-1)*2+2*(m-1)+1,:)=[A C B];
		index((n-1)*(ny-1)*2+2*m,:)=[B C D];
	end
end

%Scale  x and y
x=x/max(x)*Lx;
y=y/max(y)*Ly;

%create segments
segments=zeros(2*(nx-1)+2*(ny-1),3);
%left edge:
segments(1:ny-1,:)=[[2:ny]' [1:ny-1]' 2*[1:ny-1]'-1];
%right edge:
segments(ny:2*(ny-1),:)=[[ny*(nx-1)+1:nx*ny-1]' [ny*(nx-1)+2:nx*ny]' 2*[(ny-1)*(nx-2)+1:(nx-1)*(ny-1)]'];
%front edge:
segments(2*(ny-1)+1:2*(ny-1)+(nx-1),:)=[[2*ny:ny:ny*nx]' [ny:ny:ny*(nx-1)]' [2*(ny-1):2*(ny-1):2*(nx-1)*(ny-1)]'];
%back edge
segments(2*(ny-1)+(nx-1)+1:2*(nx-1)+2*(ny-1),:)=[[1:ny:(nx-2)*ny+1]' [ny+1:ny:ny*(nx-1)+1]' [1:2*(ny-1):2*(nx-2)*(ny-1)+1]'];

%plug coordinates and nodes
md.mesh=mesh2d();
md.mesh.x=x;
md.mesh.y=y;
md.mesh.numberofvertices=nods;
md.mesh.vertexonboundary=zeros(nods,1);md.mesh.vertexonboundary(segments(:,1:2))=1;

%plug elements
md.mesh.elements=index;
md.mesh.segments=segments;
md.mesh.numberofelements=nel;

%Now, build the connectivity tables for this mesh.
md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
