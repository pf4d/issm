function applyoptions(md, data, datatype, options, canvas, gl, node){ //{{{
	//APPLYOPTIONS - apply colobar, text, cloud, and expdisp options to current plot
	//
	//   Usage:
	//      applyoptions(md, data, options)
	//
	//   See also: PLOTMODEL, PARSE_OPTIONS
	
	//{{{ colorbar
	if (options.exist('colorbar')) {
		if (options.getfieldvalue('colorbar')==1) {
			//{{{ Handle movie data
			if (datatype == 5) {
				data = data[0];
			} //}}}
			//{{{ Variable options initialization
			var caxis = options.getfieldvalue('caxis',[ArrayMin(data), ArrayMax(data)]);
			var colorbarinnerlabels = options.getfieldvalue('colorbarinnerlabels','off');
			var ccanvasid, ctitleid, clabelsid, ccanvas, ctitle, clabels, ccontext, cmap, colorbar, cwidth, cheight, cgradient, color, y, x;
			//}}}
			//{{{ Create colorbar labels 
			var labels = [];
			var cdivisions = options.getfieldvalue('colorbarnticks', 6);
			var caxisdelta = caxis[1] - caxis[0];
			var precision = options.getfieldvalue('colorbarprecision', 3);
			if (options.getfieldvalue('log','off')!='off') {
				for (var i=cdivisions; i >= 0; i--) {
					var scale = (Math.log10(caxis[1])-Math.log10(caxis[0]))/Math.log10(options.getfieldvalue('log', 10));
					labels[i] = (Math.pow(options.getfieldvalue('log', 10), Math.log10(caxis[0])/Math.log10(options.getfieldvalue('log', 10))+scale*(cdivisions-i)/cdivisions)).toFixed(precision);
				}
			} else {
				for (var i=cdivisions; i >= 0; i--) {
					labels[i] = (caxisdelta*(cdivisions-i)/cdivisions+caxis[0]).toFixed(precision);
				}
			} //}}}
			//{{{ Initialize colorbar canvas
			ccanvasid = options.getfieldvalue('colorbarid', options.getfieldvalue('canvasid').replace('canvas','colorbar-canvas'));			
			ccanvas = $('#'+ccanvasid)[0];
			cwidth = ccanvas.width*options.getfieldvalue('colorbarwidth', 1);
			cheight = ccanvas.height*options.getfieldvalue('colorbarheight', 1);
			ccontext = ccanvas.getContext('2d');
			ccontext.clearRect(0, 0, cwidth, cheight);
			ccontext.beginPath();
			cmap = options.getfieldvalue('colormap','jet');
			colorbar = colorbars[cmap];
			cgradient = ccontext.createLinearGradient(0, 0, 0, cheight);
			//}}}
			//{{{ Draw colorbar gradient
			for (var i=0; i < colorbar.length; i++) {
				color = colorbar[colorbar.length-i-1];
				color = [Math.round(color[0]*255), Math.round(color[1]*255), Math.round(color[2]*255)];	
				cgradient.addColorStop(i/colorbar.length*(cdivisions/(cdivisions+1.0))+(1.0/(cdivisions+1.0)),'rgba('+color.toString()+', 1.0)');
			}
			ccontext.fillStyle=cgradient;
			ccontext.fillRect(0, 0, cwidth, cheight);
			//}}}
			//{{{ Draw colorbar border
			ccontext.beginPath();
			ccontext.lineWidth='1';
			ccontext.strokeStyle=options.getfieldvalue('colorbarfontcolor','black');
			ccontext.rect(0, 0, cwidth, cheight);
			ccontext.stroke();
			//}}}
			//{{{ Draw colorbar labels
			clabelsid = options.getfieldvalue('colorbarid', ccanvasid).replace('canvas','labels');
			clabels = $('#'+clabelsid);
			if (colorbarinnerlabels=='on') {
				clabels.removeClass('sim-colorbar-labels-outer');
				clabels.addClass('sim-colorbar-labels-inner');
			}
			else {
				clabels.removeClass('sim-colorbar-labels-inner');
				clabels.addClass('sim-colorbar-labels-outer');
			}
			var clabelstring = '';
			clabels.empty();
			for (var i=0; i <= cdivisions; i++) {
				y = (i+0.5)/(cdivisions+1)*cheight;
				x = 0.2*cwidth;
				clabelstring += '<li><span>'+labels[i]+'</span></li>';
				ccontext.beginPath();
				ccontext.moveTo(0, y);
				ccontext.lineTo(x, y);
				ccontext.moveTo(cwidth-x, y);
				ccontext.lineTo(cwidth, y);
				ccontext.stroke();
			}
			clabels.append(clabelstring);
			//}}}
			//{{{ Draw colorbar title
			ctitleid = options.getfieldvalue('colorbarid', ccanvasid).replace('canvas','heading');
			ctitle = $('#'+ctitleid);
			if (options.exist('colorbartitle')) { ctitle.html(options.getfieldvalue('colorbartitle')); }
			//}}}
		} 
	} //}}}
	//{{{ texture canvas
	var tcontext, tcanvas, tcanvasid, tURL, tgradient;
	tcanvasid = 'texturecanvas';
	var tcanvas = document.getElementById(tcanvasid);
	if (tcanvas == null) {
		$('<canvas id="texturecanvas" width="256" height="256" style="display: none;"></canvas>').insertAfter('#'+String(options.getfieldvalue('canvasid')));
		tcanvas = document.getElementById(tcanvasid);
	}
	tcontext = tcanvas.getContext('2d');
	tgradient = tcontext.createLinearGradient(0, 0, 0, 256);
		
	var cmap = options.getfieldvalue('colormap','jet');
	var colorbar = colorbars[cmap];
	for (var i=0; i < colorbar.length; i++) {
		color = colorbar[colorbar.length-i-1];
		color = [Math.round(color[0]*255), Math.round(color[1]*255), Math.round(color[2]*255)];	
		tgradient.addColorStop(i/colorbar.length,'rgba('+color.toString()+', 1.0)');
	}
	
	tcontext.fillStyle = tgradient;
	tcontext.fillRect(0, 0, 256, 256);
	tURL = tcanvas.toDataURL();
	node.texture = initTexture(gl, tURL);
	node.textureCanvas = tcanvas;
	node.caxis = options.getfieldvalue('caxis',[ArrayMin(data), ArrayMax(data)]);
	//}}}
	//{{{ expdisp contours
	if (options.exist('expdisp')) {
		canvas.nodes.expdisp = Node(gl, options);
		var node = canvas.nodes.expdisp;
		
		//declare variables:  {{{
		var vertices = [];
		var indices = [];
		var colors = [];
		var rgbcolor = [];
		var xmin, xmax;
		var ymin, ymax;
		var zmin, zmax;
		var scale;
		
		//Process data and model
		var x = options.getfieldvalue('expdisp').x;
		var y = options.getfieldvalue('expdisp').y;
		var z = Array.apply(null, Array(x.length)).map(Number.prototype.valueOf, 0);
		
		if (options.getfieldvalue('expdisp').z) {
			z = options.getfieldvalue('expdisp').z;
		}
		//}}}

		//Compute coordinates and data range: //{{{
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
		//}}}

		//Compute scaling: //{{{
		var scale = 1 / (xmax - xmin);
		node.shaderName = 'colored';
		node.shader = gl.shaders[node.shaderName].program;
		node.scale = [scale, scale, scale*options.getfieldvalue('heightscale', 1)];
		node.translation = [(xmin + xmax) / (-2 / scale), (ymin + ymax) / (-2 / scale), (zmin + zmax) / (-2 / scale)];
		node.modelMatrix = updateModelMatrix(node);
		node.drawMode = gl.LINE_LOOP;
		node.drawOrder = 0;
		node.useIndexBuffer = false;
		node.disableDepthTest = true;
		//}}}

		//some defaults:
		colors.itemSize = 4;

		//retrieve some options
		var linewidth=options.getfieldvalue('linewidth', 1);
		var edgecolor=options.getfieldvalue('edgecolor','black'); //RGBCOLOR?

		vertices.itemSize = 3;
		for(var i=0; i < x.length; i++){
			vertices[vertices.length] = x[i];
			vertices[vertices.length] = y[i];
			vertices[vertices.length] = z[i];

			//edgecolor
			rgbcolor = [0.0, 0.0, 0.0];
			colors[colors.length] = rgbcolor[0];
			colors[colors.length] = rgbcolor[1];
			colors[colors.length] = rgbcolor[2];
			colors[colors.length] = 1.0;
		}

		//Initalize buffers:
		node.arrays = [vertices, colors];
		node.buffers = initBuffers(gl, node.arrays);
	} //}}}
	//{{{ cloud of points
	if (options.exist('cloud')) {
		canvas.nodes.cloud = Node(gl, options);
		var node = canvas.nodes.cloud;

		//declare variables:  {{{
		var vertices = [];
		var indices = [];
		var colors = [];
		var rgbcolor = [];
		var xmin, xmax;
		var ymin, ymax;
		var zmin, zmax;
		var scale;
		
		//Process data and model
		var x = options.getfieldvalue('cloud').x;
		var y = options.getfieldvalue('cloud').y;
		var z = Array.apply(null, Array(x.length)).map(Number.prototype.valueOf, 0);
		
		if (options.getfieldvalue('cloud').z) {
			z = options.getfieldvalue('cloud').z;
		}
		//}}}

		//Compute coordinates and data range: //{{{
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
		//}}}

		//Compute scaling: //{{{
		var scale = 1 / (xmax - xmin);
		node.shaderName = 'colored';
		node.shader = gl.shaders[node.shaderName].program;
		node.scale = [scale, scale, scale*options.getfieldvalue('heightscale', 1)];
		node.translation = [(xmin + xmax) / (-2 / scale), (ymin + ymax) / (-2 / scale), (zmin + zmax) / (-2 / scale)];
		node.modelMatrix = updateModelMatrix(node);
		node.drawMode = gl.POINTS;
		node.drawOrder = 0;
		node.useIndexBuffer = false;
		node.disableDepthTest = true;
		//}}}

		//some defaults:
		colors.itemSize = 4;

		//retrieve some options
		var linewidth=options.getfieldvalue('linewidth', 1);
		var edgecolor=options.getfieldvalue('edgecolor','black'); //RGBCOLOR?

		vertices.itemSize = 3;
		for(var i=0; i < x.length; i++){
			vertices[vertices.length] = x[i];
			vertices[vertices.length] = y[i];
			vertices[vertices.length] = z[i];

			//edgecolor
			rgbcolor = [0.0, 0.0, 0.0];
			colors[colors.length] = rgbcolor[0];
			colors[colors.length] = rgbcolor[1];
			colors[colors.length] = rgbcolor[2];
			colors[colors.length] = 1.0;
		}

		//Initalize buffers:
		node.arrays = [vertices, colors];
		node.buffers = initBuffers(gl, node.arrays);
	} //}}}
	//{{{ text display
	if (options.exist('textlabels')) {
		var textcanvas, textcanvasid;	
		textcanvasid = options.getfieldvalue('textcanvasid', options.getfieldvalue('canvasid')+'-text');
		textcanvas = $('#'+textcanvasid);
		textcanvas.textlabels = options.getfieldvalue('textlabels',[]);
		
		//setup drawing function for text canvas draw calls
		textcanvas.draw = function(canvas) {
			var textcontext, textlabels, textlabel, textcanvaswidth, textcanvasheight, textcoordinates;	
			var textposition = vec3.create();
			var mvpMatrix = mat4.create();
			
			//ensure correct canvas coordinate scaling
			textcanvaswidth = this[0].clientWidth;
			textcanvasheight = this[0].clientHeight;
			this[0].width  = textcanvaswidth;
			this[0].height = textcanvasheight;
			
			textcontext = this[0].getContext('2d');
			textlabels = options.getfieldvalue('textlabels',[]);
			textcontext.clearRect(0, 0, textcanvaswidth, textcanvasheight);
			
			//worldspace to screenspace transformation for text
			for (text in textlabels) {
				textlabel = textlabels[text];
				mat4.multiply(mvpMatrix, canvas.camera.vpMatrix, canvas.nodes.overlay.modelMatrix);
				textposition = vec3.transformMat4(textposition, textlabel.pos, mvpMatrix);
				if (textposition[2] > 1) { //clip coordinates with z > 1
					continue;
				}
				textcoordinates = [(textposition[0]+1.0)/2.0*textcanvaswidth, (-textposition[1]+1.0)/2.0*textcanvasheight]; //NDC to screenspace
				textcontext.font = String(options.getfieldvalue('colorbarfontsize', 18))+'px "Lato", Helvetica, Arial, sans-serif';
				textcontext.fillStyle = options.getfieldvalue('colorbarfontcolor','black');
				textcontext.strokeStyle = options.getfieldvalue('colorbarfontcolor','black');
				textcontext.textAlign = 'center';
				textcontext.textBaseline = 'middle';
				textcontext.fillText(textlabel.text, textcoordinates[0], textcoordinates[1]);
				textcontext.strokeText(textlabel.text, textcoordinates[0], textcoordinates[1]);
			}
		}
		canvas.textcanvas = textcanvas;
	} //}}}
	//{{{ lat long overlay
	if (options.exist('latlongoverlay')) {
		var overlaycanvasid = options.getfieldvalue('latlongoverlayid', options.getfieldvalue('canvasid')+'-overlay');
		var overlaycanvas = $('#'+overlaycanvasid)[0];
		var latitudes = {
			//"-90": 1,
			//"-65": .999,
			"-60": 0.994046875,
			//"-55": 0.983187500000002,
			//"-50": 0.97173550854167,
			"-45": 0.955729166666666,
			//"-40": 0.94218750000000218,
			//"-35": 0.94218750000000218,
			"-30": 0.9226562500000024,
			//"-25": 0.87934895833333526,
			//"-20": 0.856572916666669,
			//"-15": 0.830729166666665,
			//"-10": 0.803552708333336,
			//"-5": 0.77395833333333541,
			"0": 0.74218749999999811,
			//"5": 0.70950364583333347,
			//"10": 0.67479166666666823,
			//"15": 0.63932291666666663,
			//"20": 0.60171875,
			//"25": 0.563453125,
			"30": 0.523390625000001,
			//"35": 0.48401875,
			//"40": 0.44296875,
			"45": 0.4020001,
			//"50": 0.3578125,
			//"55": 0.311875,
			"60": 0.26953124999999978,
			//"65": 0.225390625,
			//"70": 0.18125,
			//"75": 0.13541666666666671,
			//"80": 0.08953125,
			//"85": 0.046250000000000013,
			//"90": 0.0,
		}
		var longitudes = [-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180];
		overlaycanvas.draw = function(canvas) {
			var rect = overlaycanvas.getBoundingClientRect();
			overlaycanvas.width  = rect.width;
			overlaycanvas.height = rect.height;
			var ctx = overlaycanvas.getContext('2d');
			var centerx = overlaycanvas.width / 2;
			var centery = overlaycanvas.height / 2;
			var radius = (overlaycanvas.height) / 2;
			ctx.setLineDash([5, 10]);
			for(latitude in latitudes) {
				ctx.beginPath();
				ctx.arc(centerx, centery, radius * latitudes[latitude], 0, 2 * Math.PI);
				ctx.stroke();
				ctx.font = String(options.getfieldvalue('colorbarfontsize', 18))+'px "Lato", Helvetica, Arial, sans-serif';
				ctx.fillStyle = options.getfieldvalue('colorbarfontcolor','black');
				ctx.strokeStyle = options.getfieldvalue('colorbarfontcolor','black');
				ctx.textAlign = 'center';
				ctx.textBaseline = 'middle';
				ctx.fillText(latitude, centerx, centery + radius * latitudes[latitude]);
				ctx.strokeText(latitude, centerx, centery + radius * latitudes[latitude]);
			}
			ctx.setLineDash([1, 0]);
			for (longitude in longitudes) {
				ctx.beginPath();
				ctx.moveTo(centerx, centery);
				ctx.lineTo(centerx + radius * Math.sin(longitudes[longitude] * DEG2RAD), centery + radius * Math.cos(longitudes[longitude] * DEG2RAD));
				ctx.stroke();
			}
		}
		canvas.overlaycanvas = overlaycanvas;
	} //}}}
	//{{{ additional rendering nodes
	if (options.exist('render')) {
		var meshresults = processmesh(md, data, options);
		var x = meshresults[0]; 
		var y = meshresults[1]; 
		var z = meshresults[2]; 
		var elements = meshresults[3];
		var is2d = meshresults[4]; 
		var isplanet = meshresults[5];
		
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
		
		var global = vec3.length([(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2]) < 6371000/10; //tolerance for global models = center is 637100 meters away from center of earth
		var atmosphereScale = 1.25;
		var translation = global ? [(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2] : [(xmin + xmax) / 2, (ymin + ymax) - 6371000, (zmin + zmax) / 2];
		
		if (options.getfieldvalue('render',[]).indexOf('sky')!=-1) {	
			//atmosphere
			var node = Node(gl);
			node.name = "atmosphere";
			node.shaderName = "SkyFromSpace";
			node.shader = gl.shaders[node.shaderName];
			node.drawOrder = 1;
			node.cullFace = gl.FRONT;
			node.enableCullFace = true;
			node.mesh = GL.Mesh.icosahedron({size: 6371000*atmosphereScale, subdivisions: 6});
			node.rotation = [0, 0, 0];
			node.translation = translation;
			node.center = [0, 0, 0];
			updateModelMatrix(node);
			canvas.nodes[node.name] = node;
		}
		if (options.getfieldvalue('render',[]).indexOf('space')!=-1) {	
			//skysphere
			node = Node(gl);
			node.name = "skysphere";
			node.shaderName = "Textured";
			node.shader = gl.shaders[node.shaderName];
			node.drawOrder = 2;
			node.cullFace = gl.FRONT;
			node.enableCullFace = true;
			node.mesh = GL.Mesh.sphere({size: 6371000*20});
			node.texture = initTexture(gl, canvas.rootPath+'textures/TychoSkymapII_t4_2k.jpg');
			node.rotation = [0, 0, 0];
			node.translation = translation;
			node.center = [0, 0, 0];
			updateModelMatrix(node);
			canvas.nodes[node.name] = node;
		}
		if (canvas.clouds.enabled) {
			//clouds
			for (var i = 0; i < canvas.clouds.quantity; i++) {
				node = Node(gl);
				node.name = "clouds" + i;
				node.shaderName = "Clouds";
				node.shader = gl.shaders[node.shaderName];
				node.drawOrder = 2;
				node.cullFace = gl.BACK;
				node.enableCullFace = true;
				node.mesh = GL.Mesh.fromURL(canvas.rootPath+'obj/cloud.obj');
				node.rotation = [0, 0, 0];
				node.scale = [2500, 2500, 2500];
				node.translation = [translation[0], translation[1] - 405000, translation[2]];
				node.center = [0, 0, 0];
				node.animation = {"time": Date.now(),"target": node.translation,"current": node.translation};
				updateModelMatrix(node);
				canvas.nodes[node.name] = node;
				//canvas.clouds.list
			}
			//TODO: Steven, please add <canvas.clouds.quantity> total cloud nodes, randomly spread over the mesh, giving each one a new name and adding them to the canvas.clouds.list so that we can track them later.
			
		}
		if (options.getfieldvalue('render',[]).indexOf('latlong')!=-1) {	
			//latlong
			node = Node(gl);
			node.name = "clouds";
			node.shaderName = "Clouds";
			node.shader = gl.shaders[node.shaderName];
			node.drawOrder = 2;
			node.cullFace = gl.BACK;
			node.enableCullFace = true;
			node.mesh = GL.Mesh.fromURL(canvas.rootPath+'obj/cloud.obj');
			node.rotation = [0, 0, 0];
			node.scale = [2500, 2500, 2500];
			node.translation = [translation[0], translation[1] - 405000, translation[2]];
			node.center = [0, 0, 0];
			node.animation = {"time": Date.now(),"target": node.translation,"current": node.translation};
			updateModelMatrix(node);
			canvas.nodes[node.name] = node;
		}
	} //}}}
} //}}}
