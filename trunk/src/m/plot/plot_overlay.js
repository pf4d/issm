function plot_overlay(md, data, options, canvas){ //{{{
	//PLOT_OVERLAY - Function for plotting a georeferenced image.  
	//
	//   Usage:
	//      plot_overlay(md, data, options, canvas);
	//
	//   See also: PLOTMODEL, PLOT_MANAGER

	//{{{ declare variables:
	var vertices = [];
	var indices = [];
	var texcoords = [];
	var nanindices = {};
	var xmin, xmax;
	var ymin, ymax;
	var zmin, zmax;
	var matrixscale, vertexscale;

	//Process data and model
	var meshresults = processmesh(md, data, options);
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
	var node = Node(gl);
	canvas.nodes[canvas.nodes.length] = node;
	node.name = "overlay";
	node.shaderName = (options.getfieldvalue('render',[]).indexOf('ground')!=-1) ? "GroundFromSpace" : "Textured";
	node.shader = gl.shaders[node.shaderName];
	node.scale = [1, 1, matrixscale];
	node.rotation = [-90, 0, 0];
	node.translation = [0, 0, 0];
	node.center = [(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2];
	node.texture = initTexture(gl, options.getfieldvalue('overlay_image'));
	node.alpha = options.getfieldvalue('outeralpha', 1.0);
	node.drawOrder = 1;
	node.maskEnabled = options.getfieldvalue('outermask','off') == 'on';
	node.maskHeight = options.getfieldvalue('outermaskheight', 150.0);
	node.maskColor = options.getfieldvalue('outermaskcolor',[0.0, 0.0, 1.0, 1.0]);
	updateModelMatrix(node);
	
	//Handle outer radaroverlay
	if (md.radaroverlay.outerx) {
		var newelements = [];
		for (var i = 0; i < md.radaroverlay.outerindex.length; i++) {
			newelements[newelements.length] = [md.radaroverlay.outerindex[i][0] + x.length, md.radaroverlay.outerindex[i][1] + y.length, md.radaroverlay.outerindex[i][2] + z.length];
		}
		x = [].concat(x, md.radaroverlay.outerx);
		y = [].concat(y, md.radaroverlay.outery);
		z = [].concat(z, md.radaroverlay.outerheight);
		elements = [].concat(elements, newelements);
		
		//Reclaculate bounds based on otuer radaroverlay
		modelxlim = [ArrayMin(x), ArrayMax(x)];
		modelylim = [ArrayMin(y), ArrayMax(y)];
		modelzlim = [ArrayMin(z), ArrayMax(z)];
		xmin = xlim[0];
		xmax = xlim[1];
		ymin = ylim[0];
		ymax = ylim[1];
		zmin = zlim[0];
		zmax = zlim[1];
		
		node.center = [node.center[0], node.center[1], -zmax];
	}
	
	var xrange = modelxlim[1] - modelxlim[0];
	var yrange = modelylim[1] - modelylim[0];
	
	var xyz = vec3.create();
	var direction = vec3.create();
	var vertex = vec3.create();
	var magnitude;

	//generate mesh:
	for(var i = 0; i < x.length; i++){
		//Check for NaN values and remove from indices array as necessary, but preserve vertex array spacing
		if (isNaN(x[i]) || isNaN(y[i]) || isNaN(z[i])) {
			nanindices[i] = i;
			vertices[vertices.length] = vertex[0];
			vertices[vertices.length] = vertex[1];
			vertices[vertices.length] = vertex[2];
			
			texcoords[texcoords.length] = 0.0;
			texcoords[texcoords.length] = 0.0;
			continue;
		}

		if (md.mesh.classname() == 'mesh3dsurface') {
			//Scale vertices
			xyz = vec3.fromValues(x[i], y[i], z[i]);
			magnitude = vec3.length(xyz) + md.geometry.surface[i] * vertexscale;
			vec3.normalize(direction, xyz);
			vec3.scale(vertex, direction, magnitude);
			vertices[vertices.length] = vertex[0];
			vertices[vertices.length] = vertex[1];
			vertices[vertices.length] = vertex[2];
			
			texcoords[texcoords.length] = Math.atan2(vertex[1], vertex[0]) / (2 * Math.PI) + 0.5;
			texcoords[texcoords.length] = Math.asin(vertex[2] / magnitude) / Math.PI + 0.5;
		}
		else {
			//Scale vertices
			xyz = vec3.fromValues(x[i], y[i], z[i]);
			magnitude = vec3.length(xyz);
			vec3.normalize(direction, xyz);
			vec3.scale(vertex, direction, magnitude);
			vertices[vertices.length] = vertex[0];
			vertices[vertices.length] = vertex[1];
			vertices[vertices.length] = vertex[2];
			
			texcoords[texcoords.length] = (x[i] - modelxlim[0]) / xrange;
			texcoords[texcoords.length] = (y[i] - modelylim[0]) / yrange;
		}
	}
	//linearize the elements array:
	var element;
	for(var i = 0; i < elements.length; i++){
		element = [elements[i][0] - 1, elements[i][1] - 1, elements[i][2] - 1];
		if (element[0] in nanindices || element[1] in nanindices || element[2] in nanindices) continue;
		indices[indices.length] = element[0];
		indices[indices.length] = element[1];
		indices[indices.length] = element[2];
	}
	node.mesh = GL.Mesh.load({vertices: vertices, coords: texcoords, triangles: indices}, null, null, gl);
} //}}}
