function mesh=gmshplanet(varargin)
%GMSHPLANET - mesh generation for a sphere. Very specific code for gmsh. From demo/sphere.geo
%
%   Available options (for more details see ISSM website http://issm.jpl.nasa.gov/):
%
%   - radius:             radius of the planet in km
%   - resolution:         resolution in km
%   - refine:             provide mesh
%   - refinemetric:       mesh quantity to specify resolution
%
%   Returns 'mesh3dsurface' type mesh
%
%   Examples:
%      md.mesh=gmshplanet('radius',6000,'resolution',100);
%      md.mesh=gmshplanet('radius',6000,'resolution',100);

	%process options
	options=pairoptions(varargin{:});
	options=deleteduplicates(options,1);

	%recover parameters:
	radius=getfieldvalue(options,'radius')*1000;
	resolution=getfieldvalue(options,'resolution')*1000;

	%initialize mesh: 
	mesh=mesh3dsurface;

	%create .geo file:  {{{
	fid=fopen('sphere.geo','w');

	fprintf(fid,'Mesh.Algorithm = 1;\n');
	if  exist(options,'refine'),
		fprintf(fid,'Mesh.Algorithm = 7;\n');
		fprintf(fid,'Mesh.CharacteristicLengthFromPoints= 0;\n');
		fprintf(fid,'Mesh.SmoothRatio= 3;\n');
		fprintf(fid,'Mesh.RemeshAlgorithm= 1;\n');
	end
	fprintf(fid,'resolution=%g;\n',resolution);
	fprintf(fid,'radius=%g;\n',radius);
	fprintf(fid,'Point(1) = {0.0,0.0,0.0,resolution};\n');
	fprintf(fid,'Point(2) = {radius,0.0,0.0,resolution};\n');
	fprintf(fid,'Point(3) = {0,radius,0.0,resolution};\n');
	fprintf(fid,'Circle(1) = {2,1,3};\n');
	fprintf(fid,'Point(4) = {-radius,0,0.0,resolution};\n');
	fprintf(fid,'Point(5) = {0,-radius,0.0,resolution};\n');
	fprintf(fid,'Circle(2) = {3,1,4};\n');
	fprintf(fid,'Circle(3) = {4,1,5};\n');
	fprintf(fid,'Circle(4) = {5,1,2};\n');
	fprintf(fid,'Point(6) = {0,0,-radius,resolution};\n');
	fprintf(fid,'Point(7) = {0,0,radius,resolution};\n');
	fprintf(fid,'Circle(5) = {3,1,6};\n');
	fprintf(fid,'Circle(6) = {6,1,5};\n');
	fprintf(fid,'Circle(7) = {5,1,7};\n');
	fprintf(fid,'Circle(8) = {7,1,3};\n');
	fprintf(fid,'Circle(9) = {2,1,7};\n');
	fprintf(fid,'Circle(10) = {7,1,4};\n');
	fprintf(fid,'Circle(11) = {4,1,6};\n');
	fprintf(fid,'Circle(12) = {6,1,2};\n');
	fprintf(fid,'Line Loop(13) = {2,8,-10};\n');
	fprintf(fid,'Ruled Surface(14) = {13};\n');
	fprintf(fid,'Line Loop(15) = {10,3,7};\n');
	fprintf(fid,'Ruled Surface(16) = {15};\n');
	fprintf(fid,'Line Loop(17) = {-8,-9,1};\n');
	fprintf(fid,'Ruled Surface(18) = {17};\n');
	fprintf(fid,'Line Loop(19) = {-11,-2,5};\n');
	fprintf(fid,'Ruled Surface(20) = {19};\n');
	fprintf(fid,'Line Loop(21) = {-5,-12,-1};\n');
	fprintf(fid,'Ruled Surface(22) = {21};\n');
	fprintf(fid,'Line Loop(23) = {-3,11,6};\n');
	fprintf(fid,'Ruled Surface(24) = {23};\n');
	fprintf(fid,'Line Loop(25) = {-7,4,9};\n');
	fprintf(fid,'Ruled Surface(26) = {25};\n');
	fprintf(fid,'Line Loop(27) = {-4,12,-6};\n');
	fprintf(fid,'Ruled Surface(28) = {27};\n');
	fprintf(fid,'Surface Loop(29) = {28,26,16,14,20,24,22,18};\n');
	fprintf(fid,'Volume(30) = {29};\n');
	fprintf(fid,'Physical Surface(1) = {28,26,16,14,20,24,22,18};\n');
	fprintf(fid,'Physical Volume(2) = 30;\n');
	fclose(fid);
	%}}}

	if  exist(options,'refine'),
		meshini=getfieldvalue(options,'refine');
		metric=getfieldvalue(options,'refinemetric');

		%create .pos file with existing mesh and refining metric:  {{{
		fid=fopen('sphere.pos','w');

		fprintf(fid,'View "background mesh" {\n');
		for i=1:meshini.numberofelements,
			fprintf(fid,'ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n',...
			meshini.x(meshini.elements(i,1)), meshini.y(meshini.elements(i,1)), meshini.z(meshini.elements(i,1)),...
			meshini.x(meshini.elements(i,2)), meshini.y(meshini.elements(i,2)), meshini.z(meshini.elements(i,2)),...
			meshini.x(meshini.elements(i,3)), meshini.y(meshini.elements(i,3)), meshini.z(meshini.elements(i,3)),...
			metric(meshini.elements(i,1)), metric(meshini.elements(i,2)), metric(meshini.elements(i,3))...
			);
		end
		fprintf(fid,'};\n');
		
		fclose(fid);
		% }}}
	end

	%call gmsh
	if  exist(options,'refine'),
		eval(['!gmsh -tol 1e-8 -2 sphere.geo -bgm sphere.pos']);
	else
		%call gmsh
		eval(['!gmsh -tol 1e-8 -2 sphere.geo']);
	end

	%import mesh:  {{{
	fid=fopen('sphere.msh','r');

	%Get Mesh format
	A=fscanf(fid,'%s',1);
	if ~strcmp(A,'$MeshFormat'),
		error(['Expecting $MeshFormat (' A ')']);
	end

	A=fscanf(fid,'%f %i %i',[1 3]);
	A=fscanf(fid,'%s',1);
	if ~strcmp(A,'$EndMeshFormat'),
		error(['Expecting $EndMeshFormat (' A ')']);
	end

	%Nodes
	A=fscanf(fid,'%s',1);
	if ~strcmp(A,'$Nodes'),
		error(['Expecting $Nodes (' A ')']);
	end

	mesh.numberofvertices=fscanf(fid,'%i',1);
	A=fscanf(fid,'%i %f %f %f',[4 mesh.numberofvertices]);
	mesh.x = A(2,:)';
	mesh.y = A(3,:)';
	mesh.z = A(4,:)';

	A=fscanf(fid,'%s',1);
	if ~strcmp(A,'$EndNodes'),
		error(['Expecting $EndNodes (' A ')']);
	end

	%Elements
	A=fscanf(fid,'%s',1);
	if ~strcmp(A,'$Elements'),
		error(['Expecting $Elements (' A ')']);
	end
	mesh.numberofelements=fscanf(fid,'%i',1);
	A=fscanf(fid,'%i %i %i %i %i %i %i %i',[8 mesh.numberofelements]);
	mesh.elements=A(6:8,:)'; 
	A=fscanf(fid,'%s',1);
	if ~strcmp(A,'$EndElements'),
		error(['Expecting $EndElements (' A ')']);
	end
	fclose(fid); 
	%}}}

	%figure out other fields in mesh3dsurface: 
	mesh.r=sqrt(mesh.x.^2+mesh.y.^2+mesh.z.^2);
	mesh.lat = asin(mesh.z./mesh.r)/pi*180;
	mesh.long = atan2(mesh.y,mesh.x)/pi*180;

	%erase files: 
	eval(['!rm -rf sphere.geo sphere.msh sphere.pos']);

	%return mesh: 
	return;
