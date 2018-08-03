function plot_mesh(md, options, canvas) { //{{{
	//PLOT_MESH - Function for plotting wireframe mesh.
	//
	//   Usage:
	//      plot_mesh(md, options, canvas);
	//
	//   See also: PLOTMODEL, PLOT_MANAGER

	//{{{ declare variables:
	var vertices = [];
	var indices = [];
	var colors = [];
	var nanindices = {};
	var xmin, xmax;
	var ymin, ymax;
	var zmin, zmax;
	var scale, matrixscale, vertexscale;
	
	//Process data and model
	var meshresults = processmesh(md,[], options);
	var x = meshresults[0]; 
	var y = meshresults[1]; 
	var z = meshresults[2]; 
	var elements = meshresults[3]; 
	var is2d = meshresults[4]; 
	var isplanet = meshresults[5];
		
	//Compue scaling through matrices for 2d meshes and vertices for 3d meshes
	if (!md.geometry.surface) {
		md.geometry.surface=NewArrayFill(md.mesh.x.length, 0);
	}
	if (md.mesh.classname() == 'mesh3dsurface') {
		matrixscale = 1;
		vertexscale = options.getfieldvalue('heightscale', 1);
	}
	else {
		if (md.geometry.surface) {
			z=md.geometry.surface;
		}	
		matrixscale = options.getfieldvalue('heightscale', 1);
		vertexscale = 0;
	}
	//}}}

	//Compute coordinates and data range:
	var modelxlim = [ArrayMin(x), ArrayMax(x)];
	var modelylim = [ArrayMin(y), ArrayMax(y)];
	var modelzlim = [ArrayMin(z), ArrayMax(z)];
	var xlim = options.getfieldvalue('xlim', modelxlim);
	var ylim = options.getfieldvalue('ylim', modelylim);
	var zlim = options.getfieldvalue('zlim', modelzlim);
	xmin = xlim[0];
	xmax = xlim[1];
	ymin = ylim[0];
	ymax = ylim[1];
	zmin = zlim[0];
	zmax = zlim[1];

	//Compute gl variables:
	var gl = canvas.gl;
	gl.makeCurrent();
	var node = Node(gl);
	canvas.nodes[canvas.nodes.length] = node;
	node.name = "mesh";
	node.shaderName = "Colored";
	node.shader = gl.shaders[node.shaderName];
	node.lineWidth = options.getfieldvalue('linewidth', 1);
	node.scale = [1, 1, matrixscale];
	node.rotation = [-90, 0, 0];
	node.translation = [0, 0, 0];
	node.center = [(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2];
	node.drawMode = gl.LINES;
	node.drawOrder = 0;
	node.maskEnabled = options.getfieldvalue('innermask','off') == 'on';
	node.maskHeight = options.getfieldvalue('innermaskheight', 150.0)*options.getfieldvalue('heightscale', 1);
	node.maskColor = options.getfieldvalue('innermaskcolor',[0.0, 0.0, 1.0, 1.0]);
	updateModelMatrix(node);

	//retrieve some options
	var edgecolor = new RGBColor(options.getfieldvalue('edgecolor','black'));
	if (edgecolor.ok) edgecolor = [edgecolor.r/255.0, edgecolor.g/255.0, edgecolor.b/255.0, 1.0];
	else throw Error(sprintf("s%s%s\n","initWebGL error message: cound not find out edgecolor color for curent canvas ", canvas));

	//{{{ node plot
	if (elements[0].length==6){ //prisms
	}
	else if (elements[0].length==4){ //tetras
	}
	else{ //2D triangular elements
		var xyz = vec3.create();
		var direction = vec3.create();
		var vertex = vec3.create();
		var magnitude;
		var color = edgecolor;
		for(var i = 0; i < x.length; i++){
			//Check for NaN values and remove from indices array as necessary, but preserve vertex array spacing
			if (isNaN(x[i]) || isNaN(y[i]) || isNaN(z[i])) {
				nanindices[i] = i;
				vertices[vertices.length] = vertex[0];
				vertices[vertices.length] = vertex[1];
				vertices[vertices.length] = vertex[2];
				
				colors[colors.length] = color[0];
				colors[colors.length] = color[1];
				colors[colors.length] = color[2];
				colors[colors.length] = color[3];
				continue;
			}
			//Scale vertices
			xyz = vec3.fromValues(x[i], y[i], z[i]);
			magnitude = vec3.length(xyz) + md.geometry.surface[i] * vertexscale;
			vec3.normalize(direction, xyz);
			vec3.scale(vertex, direction, magnitude);
			vertices[vertices.length] = vertex[0];
			vertices[vertices.length] = vertex[1];
			vertices[vertices.length] = vertex[2];

			colors[colors.length] = color[0];
			colors[colors.length] = color[1];
			colors[colors.length] = color[2];
			colors[colors.length] = color[3];
		}
		
		//linearize the elements array: 
		var element;
		for(var i = 0; i < elements.length; i++){
			element = [elements[i][0] - 1, elements[i][1] - 1, elements[i][2] - 1];
			if (element[0] in nanindices || element[1] in nanindices || element[2] in nanindices) continue;
			indices[indices.length] = element[0];
			indices[indices.length] = element[1];
			indices[indices.length] = element[1];
			indices[indices.length] = element[2];
			indices[indices.length] = element[2];
			indices[indices.length] = element[0];
		}
	}
	//}}}
	node.mesh = GL.Mesh.load({vertices: vertices, colors: colors, triangles: indices}, null, null, gl);
} //}}}
