function plot_quiver(md, options, canvas, updateVel) { //{{{
	//PLOT_QUIVER - quiver plot with colors
	//
	//   Usage:
	//      plot_quiver(md, options, canvas)
	//
	//   See also: PLOTMODEL, PLOT_MANAGER

	//{{{ declare variables:
	var vertices = [];
	var indices = [];
	var colors = [];
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
	var v = updateVel != undefined ? updateVel.vel : md.initialization.vel;
	var vx = updateVel != undefined ? updateVel.vx : md.initialization.vx;
	var vy = updateVel != undefined ? updateVel.vy : md.initialization.vy;
		
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
	canvas.nodes["velocity"] = node;
	node.name = "quiver";
	node.shaderName = "Colored";
	node.shader = gl.shaders[node.shaderName];
	node.lineWidth = options.getfieldvalue('linewidth', 1);
	node.scale = [1, 1, matrixscale];
	node.rotation = [-90, 0, 0];
	node.translation = [0, 0, 0];
	node.center = [(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2];
	node.drawMode = gl.LINES;
	node.useIndexBuffer = false;
	node.drawOrder = 0;
	node.maskEnabled = options.getfieldvalue('innermask','off') == 'on';
	node.maskHeight = options.getfieldvalue('innermaskheight', 150.0)*options.getfieldvalue('heightscale', 1);
	node.maskColor = options.getfieldvalue('innermaskcolor',[0.0, 0.0, 1.0, 1.0]);
	updateModelMatrix(node);

	//retrieve some options
	var edgecolor=new RGBColor(options.getfieldvalue('edgecolor','black'));
	if (edgecolor.ok) edgecolor = [edgecolor.r/255.0, edgecolor.g/255.0, edgecolor.b/255.0, 1.0];
	else throw Error(sprintf("s%s%s\n","initWebGL error message: cound not find out edgecolor color for curent canvas ", canvas));

	//{{{ node plot
	if (elements[0].length==6){ //prisms
	}
	else if (elements[0].length==4){ //tetras
	}
	else{ //2D triangular elements
		var xyz = vec3.create();
		var xyz = vec3.create();
		var direction = vec3.create();
		var vertex = vec3.create();
		var vertexBase = vec3.create();
		var verticesArrow = [vec3.fromValues(0.0, 0.0, 0.0), vec3.fromValues(1.0, 0.0, 0.0), vec3.fromValues(0.667, -0.167, 0.0), vec3.fromValues(1.0, 0.0, 0.0), vec3.fromValues(0.667, 0.166, 0.0), vec3.fromValues(1.0, 0.0, 0.0)];
		var magnitude;
		var color = edgecolor;
		var scaling = options.getfieldvalue('scaling', 1);
		var scale;
		for(var i = 0; i < x.length; i++){
			//Check for NaN values and remove from indices array as necessary, but preserve vertex array spacing
			if (isNaN(x[i]) || isNaN(y[i]) || isNaN(z[i])) continue;
			//Scale vertices
			xyz = vec3.fromValues(x[i], y[i], z[i]);
			magnitude = vec3.length(xyz) + md.geometry.surface[i] * vertexscale;
			vec3.normalize(direction, xyz);
			vec3.scale(vertex, direction, magnitude);
			vec3.copy(vertexBase, vertex);
			
			scale = scaling*v[i];
			var modelMatrix = mat4.create();
			var scaleMatrix = mat4.create();
			var rotationMatrix = mat4.create();
			mat4.scale(scaleMatrix, scaleMatrix, vec3.fromValues(scale, scale, scale));
			mat4.rotate(rotationMatrix, rotationMatrix, Math.atan2(vy[i], vx[i]), [0.0, 0.0, 1.0]);
			mat4.multiply(modelMatrix, rotationMatrix, scaleMatrix);

			var temp = vec3.fromValues(0.0, 0.0, 0.0);
			for (var j = 0; j < 6; j++){
				vec3.transformMat4(vertex, verticesArrow[j], modelMatrix);
				vec3.add(vertex, vertex, vertexBase);
				vertices[vertices.length] = vertex[0];
				vertices[vertices.length] = vertex[1];
				vertices[vertices.length] = vertex[2];
				
				colors[colors.length] = color[0];
				colors[colors.length] = color[1];
				colors[colors.length] = color[2];
				colors[colors.length] = color[3];
			}
		}
	}
	//}}}
	node.mesh = GL.Mesh.load({vertices: vertices, colors: colors}, null, null, gl);
} //}}}
