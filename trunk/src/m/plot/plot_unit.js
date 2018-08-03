function plot_unit(md, data, datatype, options, canvas) { //{{{
	//PLOT_UNIT - unit plot, display data
	//
	//   Usage:
	//      plot_unit(md, data, options, canvas);
	//
	//   See also: PLOTMODEL, PLOT_MANAGER

	//{{{ declare variables: 
	//Process data and model
	var meshresults = processmesh(md, data, options);
	var x = meshresults[0]; 
	var y = meshresults[1]; 
	var z = meshresults[2]; 
	var elements = meshresults[3];
	var is2d = meshresults[4]; 
	var isplanet = meshresults[5];
	
	var vertices = new Float32Array(x.length * 3);
	var texcoords = new Float32Array(x.length * 2);
	var indices = new Uint16Array(elements.length * 3);
	var nanindices = {};
	var xmin, xmax;
	var ymin, ymax;
	var zmin, zmax;
	var datamin, datamax, datadelta;
	var matrixscale, vertexscale;
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
	var caxis;

	//Compute gl variables:
	var gl = canvas.gl;
	var node = Node(gl);
	canvas.nodes[canvas.nodes.length] = node;
	canvas.unitNode = node;
	canvas.unitData = data;
	node.name = "unit";
	node.shaderName = "Textured";
	node.shader = gl.shaders[node.shaderName];
	node.scale = [1, 1, matrixscale];
	node.rotation = [-90, 0, 0];
	node.translation = [0, 0, 0];
	node.center = [(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2];
	node.alpha = options.getfieldvalue('alpha', 1.0);
	node.drawOrder = 1;
	node.maskEnabled = options.getfieldvalue('innermask','off') == 'on';
	node.maskHeight = options.getfieldvalue('innermaskheight', 150.0);
	node.maskColor = options.getfieldvalue('innermaskcolor',[0.0, 0.0, 1.0, 1.0]);
	node.enabled = options.getfieldvalue('nodata','off') == 'off';
	updateModelMatrix(node);

	switch(datatype){
		//{{{ element plot
		case 1:
			pos=ArrayFindNot(data, NaN); //needed for element on water
			if (elements[0].length==6){ //prisms
			}
			else if (elements[0].length==4){ //tetras
			}
			else{ //2D triangular elements
			}
			break;
		//}}}
		//{{{ node plot
		case 2:
			if (elements[0].length==6){ //prisms
			}
			else if (elements[0].length==4){ //tetras
			}
			else{ //triangular elements	
				caxis = options.getfieldvalue('caxis',[ArrayMin(data), ArrayMax(data)]);
				if (options.getfieldvalue('log','off')!='off') caxis = [Math.log10(caxis[0])/Math.log10(options.getfieldvalue('log', 10)), Math.log10(caxis[1])/Math.log10(options.getfieldvalue('log', 10))];
				datamin = caxis[0];
				datamax = caxis[1];
				datadelta = datamax - datamin;

				var xyz = vec3.create();
				var direction = vec3.create();
				var vertex = vec3.create();
				var magnitude;

				for(var i = 0, vindex = 0, tindex = 0; i < x.length; i++){
					//Check for NaN values and remove from indices array as necessary, but preserve vertex array spacing
					if (isNaN(x[i]) || isNaN(y[i]) || isNaN(z[i]) || isNaN(data[i])) {
						nanindices[i] = i;
						vertices[vindex++] = vertex[0];
						vertices[vindex++] = vertex[1];
						vertices[vindex++] = vertex[2];
						
						texcoords[tindex++] = 0.5;
						texcoords[tindex++] = 0.0;
						continue;
					}

					//Scale vertices
					xyz = vec3.fromValues(x[i], y[i], z[i]);
					magnitude = vec3.length(xyz) + md.geometry.surface[i] * vertexscale;
					vec3.normalize(direction, xyz);
					vec3.scale(vertex, direction, magnitude);
					vertices[vindex++] = vertex[0];
					vertices[vindex++] = vertex[1];
					vertices[vindex++] = vertex[2];

					texcoords[tindex++] = 0.5;
					texcoords[tindex++] = clamp((data[i] - datamin) / datadelta, 0.0, 1.0);
				}

				//linearize the elements array: 
				var element;
				for(var i = 0, iindex = 0; i < elements.length; i++){
					element = [elements[i][0] - 1, elements[i][1] - 1, elements[i][2] - 1];
					if (element[0] in nanindices || element[1] in nanindices || element[2] in nanindices) continue;
					indices[iindex++] = element[0];
					indices[iindex++] = element[1];
					indices[iindex++] = element[2];
				}
			}
			node.mesh = GL.Mesh.load({vertices: vertices, coords: texcoords, triangles: indices}, null, null, gl);
			node.mesh.octree = new GL.Octree(node.mesh);
			break;
		//}}}
		//{{{ quiver plot 
		case 3:
			if (is2d){
				//plot_quiver(x, y, data(:, 1), data(:, 2), options);
			}
			else{
				//plot_quiver3(x, y, z, data(:, 1), data(:, 2), data(:, 3), options);
			}
			break;
		//}}}
		//{{{ node transient plot
		case 5:
			if (elements[0].length==6){ //prisms
			}
			else if (elements[0].length==4){//tetras
			}
			else{ //triangular elements
				var xyz = vec3.create();
				var direction = vec3.create();
				var vertex = vec3.create();
				var magnitude;
				var timestamps = data[data.length-1];
				for(var i = 0, vindex = 0, tindex = 0; i < x.length; i++){
					//Check for NaN values and remove from indices array as necessary, but preserve vertex array spacing
					if (isNaN(x[i]) || isNaN(y[i]) || isNaN(z[i]) || isNaN(data[i][0])) {
						nanindices[i] = i;
					}
					else {
						//Scale vertices
						xyz = vec3.fromValues(x[i], y[i], z[i]);
						magnitude = vec3.length(xyz) + md.geometry.surface[i] * vertexscale;
						vec3.normalize(direction, xyz);
						vec3.scale(vertex, direction, magnitude);
					}
					vertices[vindex++] = vertex[0];
					vertices[vindex++] = vertex[1];
					vertices[vindex++] = vertex[2];
				}	
				//Transpose data to obtain column addressable data matrix
				data = data[0].map(function(col, i) { 
					return data.map(function(row) { 
						return row[i]
					})
				});
				//Prevent evaluation of datasubarray min/max if caxis exists
				if (options.exist('caxis')) caxis = options.getfieldvalue('caxis');
				else caxis = [ArrayMin(data[0]), ArrayMax(data[0].slice(0,-1))];
				if (options.getfieldvalue('log','off')!='off') caxis = [Math.log10(caxis[0])/Math.log10(options.getfieldvalue('log', 10)), Math.log10(caxis[1])/Math.log10(options.getfieldvalue('log', 10))];
				//Prepare texcoords to hold array of data values
				texcoords = [];
				for(var i = 0; i < data.length; i++){					
					datamin = caxis[0];
					datamax = caxis[1];
					datadelta = datamax - datamin;
					//Precalculate arrays for each datasubarray, trimming off timestamp value by using x.length instead of data[i].length
					texcoords[i] = new Float32Array(x.length * 2);
					for(var j = 0, index = 0; j < x.length; j++){
						texcoords[i][index++] = 0.5;
						texcoords[i][index++] = clamp((data[i][j] - datamin) / datadelta, 0.0, 1.0);
					}
				}
				
				//linearize the elements array:
				var element;
				for(var i = 0, iindex = 0; i < elements.length; i++){
					element = [elements[i][0] - 1, elements[i][1] - 1, elements[i][2] - 1];
					if (element[0] in nanindices || element[1] in nanindices || element[2] in nanindices) continue;
					indices[iindex++] = element[0];
					indices[iindex++] = element[1];
					indices[iindex++] = element[2];
				}
			
				//Initialize movie loop
				node.movieLoop = canvas.animation.loop;
				node.movieInterval = 1000 / canvas.animation.fps;
				node.movieTimestamps = timestamps;
				node.movieLength = timestamps.length;
				node.movieFrame = 0;
				canvas.dataMarkers.values = [];
				var quiverVelFrames = {};
				for(var i=0; i < md.results.length; i++){
					quiverVelFrames[Math.floor(md.results[i].time)] = md.results[i];
				}

				if (canvas.animation.handler !== 0) {
					console.log("clearing...");
					clearInterval(canvas.animation.handler)
				}
				//TODO: Move this into webgl.js
				canvas.animation.handler = setInterval(function () {
					node.movieFrame = canvas.animation.frame;
					if (canvas.animation.play && canvas.animation.increment) {
						if (node.movieFrame == node.movieLength - 1) {
							if (node.movieLoop) {
								node.movieFrame = 0;
							}
							else { 
								toggleMoviePlay(canvas);
							}
						}
						else {
							node.movieFrame = node.movieFrame + 1;
						}
						if (canvas.animation.lastFrame != canvas.animation.frame) {
							updateMarker(canvas, false);
						}
					}
					
					if (canvas.progressBar) {
						canvas.progressBar.val(node.movieFrame);
						canvas.progressBar.slider('refresh');
					}
					if (canvas.timeLabel) { canvas.timeLabel.html(node.movieTimestamps[node.movieFrame].toFixed(0) + " " + options.getfieldvalue("movietimeunit","yr")); }
					
					var buffer = node.mesh.getBuffer("coords");
					buffer.data = texcoords[node.movieFrame];
					buffer.upload(canvas.gl.DYNAMIC_DRAW);
					node.mesh.octree = new GL.Octree(node.mesh);
					node.texcoords = texcoords;
					if(options.getfieldvalue('quiver') == 'on'){
						plot_quiver(md, options, canvas, {vel: quiverVelFrames[node.movieFrame].Vel, vx: quiverVelFrames[node.movieFrame].Vx, vy: quiverVelFrames[node.movieFrame].Vy});
					}
					canvas.animation.lastFrame = canvas.animation.frame;
					canvas.animation.frame = node.movieFrame;					
				}, node.movieInterval);
				
				if (canvas.progressBar) {
					canvas.animation.frame = 0;
					canvas.progressBar.val(0);
					canvas.progressBar.attr('max', node.movieLength-1);
					canvas.progressBar.slider('refresh');
				}
				
			}
			node.mesh = GL.Mesh.load({vertices: vertices, coords: texcoords[0], triangles: indices}, null, null, gl);
			node.mesh.octree = new GL.Octree(node.mesh);
			break;
		//}}}
		default:
			throw Error(sprintf("%s%i%s\n",'case ', datatype,' not supported'));
	}
} //}}}
