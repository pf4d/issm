function  vector_average=DepthAverage(md,vector)
%DEPTHAVERAGE - computes depth average of 3d vector using the trapezoidal rule, and returns the value on 2d mesh. 
%
%   Usage:
%      vector_average=DepthAverage(md,vector);
%
%   Example:
%      vel_bar=DepthAverage(md,md.initialization.vel);

%check that the model given in input is 3d
if ~strcmp(md.mesh.elementtype(),'Penta');
	error('DepthAverage error message: the model given in input must be 3d')
end

%nods data
if (length(vector)==md.mesh.numberofvertices),
	vector_average=zeros(md.mesh.numberofvertices2d,1);
	for i=1:md.mesh.numberoflayers-1,
		vector_average=vector_average+(project2d(md,vector,i)+project2d(md,vector,i+1))/2.*(project2d(md,md.mesh.z,i+1)-project2d(md,md.mesh.z,i));
	end
	vector_average=vector_average./project2d(md,md.geometry.thickness,1);

%element data
elseif (length(vector)==md.mesh.numberofelements),
	vector_average=zeros(md.mesh.numberofelements2d,1);
	for i=1:md.mesh.numberoflayers-1,
		vector_average=vector_average+project2d(md,vector,i).*(project2d(md,md.mesh.z,i+1)-project2d(md,md.mesh.z,i));
	end
	vector_average=vector_average./project2d(md,md.geometry.thickness,1);

else
	error('vector size not supported yet');
end
