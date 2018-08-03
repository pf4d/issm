/*This is where we have all our webgl relevant functionality for the plotting routines: */
//{{{ Canvas Initialization
function initCanvas(options) {
	//Initialize open Gl for each canvas and clear any previous animation handlers, once per plotmodel call:
	canvas = document.getElementById(options.getfieldvalue('canvasid'));
	//var canvas = document.getElementById(options.getfieldvalue('canvasid'));
	if (!canvas.initialized) {
		typedArraySliceSupport();
		if (!isEmptyOrUndefined(canvas.draw) && canvas.draw.handler !== 0)	{ window.cancelAnimationFrame(canvas.draw.handler); }
		if (!isEmptyOrUndefined(canvas.animation) && canvas.animation.handler !== 0) { clearInterval(canvas.animation.handler); }
		initWebGL(canvas, options);
		initializeMarker(canvas);
		canvas.nodes = [];
		draw(canvas);
		canvas.initialized = true;
	}
	return canvas;
}
function initWebGL(canvas, options) { //{{{
	//Initialize canvas.gl on page load, reusing gl context on additional runs
	var gl;
	if (!canvas.gl) {
		gl = GL.create({canvas: canvas});
		gl.enable(gl.DEPTH_TEST); // Enable depth testing
		gl.depthFunc(gl.LEQUAL); // Near things obscure far things
		gl.enable(gl.BLEND); // Enable color blending/overlay
		gl.enable(gl.CULL_FACE); // Enable face culling
		gl.cullFace(gl.FRONT);
		gl.shaders = loadShaders(gl, options.getfieldvalue('rootpath', '../../../js/')); // Load shaders and store them in gl object
		gl.textures = {};
		
		// Add event listeners for canvas
		var displayview = options.getfieldvalue('displayview', 'off') == 'on';
		var displayzoom = options.getfieldvalue('displayzoom', 'off') == 'on';
		var mc = new Hammer.Manager(canvas);
		
		mc.add( new Hammer.Tap({event: 'singletap' }) );
		mc.add(new Hammer.Pan({threshold: 0, pointers: 0}));
		mc.add(new Hammer.Pinch({threshold: 0})).recognizeWith(mc.get('pan'));
		mc.on('singletap', function (ev) {onTap(ev, canvas);});
		mc.on('panstart panmove', function (ev) {onPan(ev, canvas, displayview);});
		mc.on('pinchstart pinchmove', function (ev) {onPinch(ev, canvas, displayview);});
		
		canvas.addEventListener('mousewheel', function (ev) {onZoom(ev, canvas, displayzoom)}, false);
		canvas.addEventListener('DOMMouseScroll', function (ev) {onZoom(ev, canvas, displayzoom)}, false);
	}
	else {
		gl = canvas.gl;
	}
	
	// Add context state variables
	canvas.gl = gl;
	canvas.textcanvas = null;
	canvas.overlaycanvas = null;
	canvas.unitNode = {};
	canvas.unitData = {};
	canvas.controlSensitivity = options.getfieldvalue('controlsensitivity', 1);
	canvas.id = options.getfieldvalue('canvasid', '.sim-canvas');
	canvas.rootPath = options.getfieldvalue('rootpath', '../../../js/');
	canvas.selector = $('#' + canvas.id);
	var backgroundcolor = new RGBColor(options.getfieldvalue('backgroundcolor', 'lightcyan'));
	if (backgroundcolor.ok) { canvas.backgroundcolor = [backgroundcolor.r/255.0, backgroundcolor.g/255.0, backgroundcolor.b/255.0, 1.0]; }
	else { throw Error(sprintf('s%s%s\n','initWebGL error message: cound not find out background color for curent canvas ', canvas)); }
	
	//Property intiialization, using values from options first, then from default values.
	var animation = options.getfieldvalue('movies', {});
	canvas.animation = {
		frame: defaultFor(animation.frame, 0),
		play: defaultFor(animation.play, true),
		increment: defaultFor(animation.increment, true),
		fps: defaultFor(animation.fps, 4),
		loop: defaultFor(animation.loop, true),
		handler: defaultFor(animation.handler, 0)
	}
	var brush = options.getfieldvalue('brush', {});
	canvas.brush = {
		enabled: defaultFor(brush.enabled, false),
		strength: defaultFor(brush.strength, 0.075),
		falloff: defaultFor(brush.falloff, 0.5),
		hit: defaultFor(brush.hit, {})
	};
	var camera = options.getfieldvalue('camera', {});
	canvas.camera = {
		position: defaultFor(camera.position, vec3.create()),
		rotation: defaultFor(camera.rotation, vec3.create()),
		near: defaultFor(camera.near, 1e3),
		far: defaultFor(camera.far, 1e10),
		fov: defaultFor(camera.fov, 45),
		vMatrix: defaultFor(camera.vMatrix, mat4.create()),
		pMatrix: defaultFor(camera.pMatrix, mat4.create()),
		vpMatrix: defaultFor(camera.vpMatrix, mat4.create()),
		vInverseMatrix: defaultFor(camera.vInverseMatrix, mat4.create()),
		pInverseMatrix: defaultFor(camera.pInverseMatrix, mat4.create()),
		vpInverseMatrix: defaultFor(camera.vpInverseMatrix, mat4.create()),
		ready: defaultFor(camera.ready, false)
	};
	var clouds = options.getfieldvalue('clouds', {});
	canvas.clouds = {
		enabled: defaultFor(clouds.enabled, false),
		height: defaultFor(clouds.height, 7500),
		quantity: defaultFor(clouds.quantity, 10)
	};
	var dataMarkers = options.getfieldvalue('datamarkers', {});
	canvas.dataMarkers = {
		enabled: defaultFor(dataMarkers.enabled, true),
		values: defaultFor(dataMarkers.values, []),
		image: defaultFor(dataMarkers.image, canvas.rootPath+'textures/data_marker.svg'),
		size: defaultFor(dataMarkers.size, [32, 32]),
		format: defaultFor(dataMarkers.format, ['X: %.2e<br>Y: %.2e<br>Z: %.2e<br>Value: %0.1f', 'x', 'y', 'z', 'value']),
		animated: defaultFor(dataMarkers.animated, false),
		labels: defaultFor(dataMarkers.labels, []),
		font: defaultFor(dataMarkers.font, ''),
		marker: defaultFor(dataMarkers.marker, document.getElementById('sim-data-marker-' + canvas.id))
	};
	var draw = options.getfieldvalue('draw', {});
	canvas.draw = {
		ready: defaultFor(draw.ready, false),
		handler: defaultFor(draw.handler, null)
	};
	var view = options.getfieldvalue('view', {});
	canvas.view = {
		position: defaultFor(view.position, [0.0, 0.0, 0.0]),
		rotation: defaultFor(view.rotation, [0, 90]),
		zoom: defaultFor(view.zoom, 1.0),
		zoomLimits: defaultFor(view.zoomLimits, [0.001, 100.0]),
		lastZoom: defaultFor(view.lastZoom, 1.0),
		azimuthLimits: defaultFor(view.azimuthLimits, [0, 360]),
		elevationLimits: defaultFor(view.elevationLimits, [-180, 180]),
		panningEnabled: defaultFor(view.panningEnabled, false),
		twod: defaultFor(view.twod, false)
	};

	//Override with parameters from URL, if any
	//TODO: Make permalinks more robust and less interdependent on UI
	if (!canvas.usedparemters) {
		function getJsonFromUrl() {
			var query = location.search.substr(1);
			var result = {};
			query.split('&').forEach(function(part) {
				var item = part.split('=');
				result[item[0]] = decodeURIComponent(item[1]);
			});
			return result;
		}
		parameters = getJsonFromUrl();
		
		if (parameters['view']) {
			canvas.view = JSON.parse(parameters['view']);
		}
		if (parameters['initial']) {
			initial = JSON.parse(parameters['initial']);
			if (!initial) {
				if (typeof SolveGlacier == 'function') {
					SolveGlacier();
				}
				if (typeof SolveSlr == 'function') {
					SolveSlr();
				}
			}
		}
		canvas.usedparemters = true;
	}
} //}}}
function generatePermalink() { //{{{
	var permalink = window.location.origin + window.location.pathname + '&view=' + JSON.stringify(canvas.view) + '&initial=' + JSON.stringify(initial);
	window.prompt('Share this simulation: ', permalink);
} //}}}
function loadShaders(gl, rootPath) { //{{{
	var shaders = {};
	shaders.Colored = new GL.Shader.fromURL(rootPath+'shaders/Colored.vsh', rootPath+'shaders/Colored.fsh', null, gl);
	shaders.Textured = new GL.Shader.fromURL(rootPath+'shaders/Textured.vsh', rootPath+'shaders/Textured.fsh', null, gl);
	shaders.SkyFromSpace = new GL.Shader.fromURL(rootPath+'shaders/SkyFromSpace.vert', rootPath+'shaders/SkyFromSpace.frag', null, gl);
	shaders.GroundFromSpace = new GL.Shader.fromURL(rootPath+'shaders/GroundFromSpace.vert', rootPath+'shaders/GroundFromSpace.frag', null, gl);
	shaders.Clouds = new GL.Shader.fromURL(rootPath+'shaders/Clouds.vert', rootPath+'shaders/Clouds.frag', null, gl);
	return shaders;
} //}}}
function initTexture(gl, imageSource) { //{{{
	//Initialize textures, or load from memory if they already exist.
	if (isEmptyOrUndefined(gl.textures[imageSource])) {
		gl.textures[imageSource] = GL.Texture.fromURL(imageSource, {minFilter: gl.LINEAR_MIPMAP_LINEAR, magFilter: gl.LINEAR}, null, gl);
	}
	return gl.textures[imageSource];
} //}}}
function Node(gl) { //{{{
	//Returns a Node object that contains default display states for webgl object. center represents pivot point of rotation.
	return {
		alpha: 1.0,
		buffers: [],
		cullFace: gl.BACK,
		disableDepthTest: false, 
		drawMode: gl.TRIANGLES,
		drawOrder: 0,
		enabled: true,
		enableCullFace: true,
		hideOcean: false,
		lineWidth: 1.0,
		maskEnabled: false,
		maskHeight: 150.0,
		maskColor: vec4.fromValues(0.0, 0.0, 1.0, 1.0),
		mesh: null,
		name: 'node',
		shaderName: 'Colored',
		shader: gl.shaders.Colored,
		texture: null,
		useIndexBuffer: true,
		center: vec3.create(), 
		scale: vec3.fromValues(1, 1, 1),
		rotation: vec3.create(),
		translation: vec3.create(),
		modelMatrix: mat4.create(),
		rotationMatrix: mat4.create(),
		inverseModelMatrix: mat4.create(),
		inverseRotationMatrix: mat4.create()
	};
} //}}}
function debugNodes(canvasid) { //{{{
	var canvasid = canvasid || '.sim-canvas';
	var nodes = document.getElementById(canvasid).nodes;
	console.log(canvasid, 'Nodes:');
	for (var node in nodes) {
		console.log('name: ', nodes[node].name, ' node: ', nodes[node], ' mesh: ', nodes[node].mesh, ' translation: ', nodes[node].translation, ' center:', nodes[node].center, ' rotation:', nodes[node].rotation);
	}
	return nodes;
} //}}}
function updateModelMatrix(node) { //{{{
	var modelMatrix = mat4.create();

	var translationMatrix = mat4.create();
	mat4.translate(translationMatrix, translationMatrix, vec3.negate(vec3.create(), node.center)); //scale/rotation centering
	mat4.multiply(modelMatrix, translationMatrix, modelMatrix);
	
	var scaleMatrix = mat4.create();
	mat4.scale(scaleMatrix, scaleMatrix, node.scale);
	mat4.multiply(modelMatrix, scaleMatrix, modelMatrix);
	
	var rotationMatrix = mat4.create();
	var zRotationMatrix = mat4.create();	
	mat4.rotate(zRotationMatrix, zRotationMatrix, DEG2RAD * node.rotation[2], [0.0, 0.0, 1.0]);
	mat4.multiply(rotationMatrix, zRotationMatrix, rotationMatrix);
	var yRotationMatrix = mat4.create();	
	mat4.rotate(yRotationMatrix, yRotationMatrix, DEG2RAD * node.rotation[1], [0.0, 1.0, 0.0]);
	mat4.multiply(rotationMatrix, yRotationMatrix, rotationMatrix);
	var xRotationMatrix = mat4.create();	
	mat4.rotate(xRotationMatrix, xRotationMatrix, DEG2RAD * node.rotation[0], [1.0, 0.0, 0.0]);
	mat4.multiply(rotationMatrix, xRotationMatrix, rotationMatrix);
	mat4.multiply(modelMatrix, rotationMatrix, modelMatrix);	
	
	mat4.identity(translationMatrix);
	mat4.translate(translationMatrix, translationMatrix, node.center); //relative translation
	mat4.multiply(modelMatrix, translationMatrix, modelMatrix);
	
	mat4.identity(translationMatrix);
	mat4.translate(translationMatrix, translationMatrix, node.translation); //absolute translation
	mat4.multiply(modelMatrix, translationMatrix, modelMatrix);
	
	node.modelMatrix = modelMatrix;
	node.inverseModelMatrix = mat4.invert(mat4.create(), modelMatrix);
	node.rotationMatrix = rotationMatrix;
	node.inverseRotationMatrix = mat4.invert(mat4.create(), rotationMatrix);
} //}}}
function clamp(value, min, max) { //{{{
	return Math.max(min, Math.min(value, max));
} //}}}
function defaultFor(name, value) { //{{{
	return typeof name !== 'undefined' ? name : value;
} //}}}
function isEmptyOrUndefined(object) { //{{{
	return object === undefined || Object.getOwnPropertyNames(object).length === 0;
} //}}}
function recover(canvasid, name, value) { //{{{
	//Traverse canvas object tree for property defined by dot delimited string, returning it, or a default value if it is not found.
	var object = document.getElementById(canvasid);
	var properties = name.split('.');
	for (var i = 0; i < properties.length; ++i) {
		object = object[properties[i]];
		if (typeof object === 'undefined') { break; }
    }
	return defaultFor(object, value);
} //}}}
function typedArraySliceSupport() { //{{{
	//TypedArray compatibility for Safari/IE
	if (typeof Int8Array !== 'undefined') {
		if (!Int8Array.prototype.fill) { Int8Array.prototype.fill = Array.prototype.fill; }
		if (!Int8Array.prototype.slice) { Int8Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Uint8Array !== 'undefined') {
		if (!Uint8Array.prototype.fill) { Uint8Array.prototype.fill = Array.prototype.fill; }
		if (!Uint8Array.prototype.slice) { Uint8Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Uint8ClampedArray !== 'undefined') {
		if (!Uint8ClampedArray.prototype.fill) { Uint8ClampedArray.prototype.fill = Array.prototype.fill; }
		if (!Uint8ClampedArray.prototype.slice) { Uint8ClampedArray.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Int16Array !== 'undefined') {
		if (!Int16Array.prototype.fill) { Int16Array.prototype.fill = Array.prototype.fill; }
		if (!Int16Array.prototype.slice) { Int16Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Uint16Array !== 'undefined') {
		if (!Uint16Array.prototype.fill) { Uint16Array.prototype.fill = Array.prototype.fill; }
		if (!Uint16Array.prototype.slice) { Uint16Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Int32Array !== 'undefined') {
		if (!Int32Array.prototype.fill) { Int32Array.prototype.fill = Array.prototype.fill; }
		if (!Int32Array.prototype.slice) { Int32Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Uint32Array !== 'undefined') {
		if (!Uint32Array.prototype.fill) { Uint32Array.prototype.fill = Array.prototype.fill; }
		if (!Uint32Array.prototype.slice) { Uint32Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Float32Array !== 'undefined') {
		if (!Float32Array.prototype.fill) { Float32Array.prototype.fill = Array.prototype.fill; }
		if (!Float32Array.prototype.slice) { Float32Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Float64Array !== 'undefined') {
		if (!Float64Array.prototype.fill) { Float64Array.prototype.fill = Array.prototype.fill; }
		if (!Float64Array.prototype.slice) { Float64Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof TypedArray !== 'undefined') {
		if (!TypedArray.prototype.fill) { TypedArray.prototype.fill = Array.prototype.fill; }
		if (!TypedArray.prototype.slice) { TypedArray.prototype.slice = Array.prototype.slice; }
	}
} //}}}
//}}}
//{{{ Interface Functions
function onTap(ev, canvas) { //{{{
	//Sets up a marker on a canvas that will track a point on the mesh. Can be dismissed by closing the display or clicking the marker.
	ev.preventDefault();
	if (!canvas.dataMarkers.enabled) { return; }
	var hit = raycast(canvas, ev.srcEvent.layerX, ev.srcEvent.layerY);
	canvas.dataMarkers.marker.hit = hit;
	canvas.brush.hit = hit;
	updateMarker(canvas, true);
	brushModify(canvas);
} //}}}
function onPan(ev, canvas, displaylog) { //{{{
	ev.preventDefault();
	
	if (canvas.dataMarkers.enabled == 'on') {
		canvas.brush.hit = raycast(canvas, ev.srcEvent.layerX, ev.srcEvent.layerY);
		brushModify(canvas);
	}
	
	if (ev.type == 'panstart') {
		canvas.lastDeltaX = 0;
		canvas.lastDeltaY = 0;
	}
	if (ev.srcEvent.shiftKey || ev.pointers.length == 2) {
		if (!canvas.view.panningEnabled) return;
		var deltaX = (canvas.lastDeltaX - ev.deltaX) / canvas.clientWidth / canvas.view.zoom * 2 * canvas.controlSensitivity * 6.371e6;
		var deltaY = (canvas.lastDeltaY - ev.deltaY) / canvas.clientHeight / canvas.view.zoom * 2 * canvas.controlSensitivity * 6.371e6;
		
		if (canvas.view.twod) {
			canvas.view.position[0] += Math.cos(DEG2RAD * canvas.view.rotation[0]) * deltaX - Math.sin(DEG2RAD * 0) * deltaY;
			canvas.view.position[2] += Math.sin(DEG2RAD * canvas.view.rotation[0]) * deltaX + Math.cos(DEG2RAD * 0) * deltaY;
		}
		else {
			canvas.view.position[0] += Math.cos(DEG2RAD * canvas.view.rotation[0]) * deltaX - Math.sin(DEG2RAD * canvas.view.rotation[0]) * deltaY;
			canvas.view.position[2] += Math.sin(DEG2RAD * canvas.view.rotation[0]) * deltaX + Math.cos(DEG2RAD * canvas.view.rotation[0]) * deltaY;
		}
	}
	
	else {
		canvas.view.rotation[0] += (canvas.lastDeltaX - ev.deltaX) / canvas.clientWidth * -2 * canvas.controlSensitivity * RAD2DEG;
		canvas.view.rotation[1] += (canvas.lastDeltaY - ev.deltaY) / canvas.clientHeight * -2 * canvas.controlSensitivity * RAD2DEG;
		
		if (canvas.view.rotation[0] > 360) { canvas.view.rotation[0] -= 360; };
		if (canvas.view.rotation[0] < -360) { canvas.view.rotation[0] += 360; };
		if (canvas.view.rotation[1] > 180) { canvas.view.rotation[1] -= 360; };
		if (canvas.view.rotation[1] < -180) { canvas.view.rotation[1] += 360; };
		
		canvas.view.rotation[0] = clamp(canvas.view.rotation[0], canvas.view.azimuthLimits[0], canvas.view.azimuthLimits[1]);
		canvas.view.rotation[1] = clamp(canvas.view.rotation[1], canvas.view.elevationLimits[0], canvas.view.elevationLimits[1])
	}
	canvas.lastDeltaX = ev.deltaX;
	canvas.lastDeltaY = ev.deltaY;
	
	repositionMarker(canvas);
	
	if (displaylog) { console.log(canvas.view.rotation); }
} //}}}
function onPinch(ev, canvas, displaylog) { //{{{
	ev.preventDefault();
	if (ev.type == 'pinchstart') { canvas.view.lastZoom = canvas.view.zoom; }
	else { modifyZoom(ev.scale * canvas.view.lastZoom, canvas, displaylog); }
} //}}}
function onZoom(ev, canvas, displaylog) { //{{{
	ev.preventDefault();
	var delta = clamp(ev.scale || ev.wheelDelta || -ev.detail, -1, 1) * canvas.controlSensitivity * canvas.view.zoom / 20;
	modifyZoom(canvas.view.zoom + delta, canvas, displaylog);
} //}}}
function modifyZoom(value, canvas, displaylog) { //{{{
	canvas.view.zoom = clamp(value, canvas.view.zoomLimits[0], canvas.view.zoomLimits[1]);
	repositionMarker(canvas);
	if (displaylog) { console.log(canvas.view.zoom); }
} //}}}
function modifyDataMarkersEnabled(value, canvas) { //{{{
	canvas.dataMarkers.enabled = value;
} //}}}
function toggleMoviePlay(canvas) { //{{{
	canvas.animation.play = !canvas.animation.play;
	if (canvas.animation.play){
		canvas.playButton.find('span').removeClass('fa-play');
		canvas.playButton.find('span').addClass('fa-pause');
	}
	else{
		canvas.playButton.find('span').removeClass('fa-pause');
		canvas.playButton.find('span').addClass('fa-play');
	}
} //}}}
function onSlideStart(canvas, progressBar) { //{{{
	if (!isEmptyOrUndefined(canvas.animation)) {
		canvas.animation.increment = false;	
		canvas.animation.frame = parseInt($(progressBar).val());
		//console.log(canvas.animation.frame);
		//updateMarker(canvas, false);
	}
} //}}}
function onSlideChange(canvas, progressBar) { //{{{
	if (!isEmptyOrUndefined(canvas.animation)) {
		canvas.animation.frame = parseInt($(progressBar).val());
		//console.log("change");
		updateMarker(canvas, false);
	}
} //}}}
function onSlideStop(canvas, progressBar) { //{{{
	if (!isEmptyOrUndefined(canvas.animation)) {
		canvas.animation.increment = true;	
		canvas.animation.frame = parseInt($(progressBar).val());
		//console.log(canvas.animation.frame);
		//updateMarker(canvas, false);
	}
} //}}}
//}}}
//{{{ Interaction Functions
function raycast(canvas, x, y) { //{{{
	//Performs raycast on canvas.unitNode.mesh using x/y screen coordinates. Returns hit objects with hit position, coords, and indicies of ray-triangle intersection.
	//TODO: Diagnose marker issues with orthographic views and slr-eustatic updates when switching between basins.
	var inverseMVPMatrix = mat4.invert(mat4.create(), mat4.multiply(mat4.create(), canvas.camera.vpMatrix, canvas.unitNode.modelMatrix));
	var origin = vec3.transformMat4(vec3.create(), [(x - canvas.width / 2) / (canvas.width / 2), (canvas.height / 2 - y) / (canvas.height / 2), 0], inverseMVPMatrix);
	var far = far || vec3.transformMat4(vec3.create(), [(x - canvas.width / 2) / (canvas.width / 2), (canvas.height / 2 - y) / (canvas.height / 2), 1.0], inverseMVPMatrix);
	var ray = vec3.subtract(vec3.create(), far, origin);

	var mesh = canvas.unitNode.mesh;
	if (!mesh) { return; }
	if (!mesh.octree) { mesh.octree = new GL.Octree(mesh); }
	
	var hit = mesh.octree.testRay(origin, ray, 1e3, 1e10);
	
	if(!hit) { return; }
	
	hit.modelPos = vec3.copy(vec3.create(), hit.pos);
	vec3.transformMat4(hit.pos, hit.pos, canvas.unitNode.modelMatrix);
	vec3.transformMat4(hit.normal, hit.normal, canvas.unitNode.modelMatrix);

	return hit;
} //}}}
function brushModify(canvas) { //{{{
	//This function takes in the canvas and x/y coordinates, performing a raycast against the mesh, and modifies the mesh using a the canvas.brush.strength and canvas.brush.falloff properties.
	//Currently the brush extends to the raycasted element and its immediate neighbors.
	//TODO: Allow variable brush size/additional neighbors. Allow the function to work on multiple models (currently hardcoded to md model).
	if (!canvas.unitNode || canvas.brush.enabled != 'on') { return; }
	
	var hit = canvas.brush.hit;

	if (hit) {
		var bufferVertices = canvas.unitNode.mesh.getBuffer('vertices');
		var vertices = bufferVertices.data;
		var bufferCoords = canvas.unitNode.mesh.getBuffer('coords');
		var coords = bufferCoords.data;

		//Query nearby elements and store indicies of affected vertices using pregenerated vertexconnectivity list (from NodeConnectivity)
		var baseIndices = new Set(hit.indices);
		var connectedIndices = new Set(hit.indices);
		var connectedElement;
		var indices;
		var lengthIndex = md.mesh.vertexconnectivity[0].length - 1;
		var length;
		for (var i = 0; i < 3; i++) {
			length = md.mesh.vertexconnectivity[hit.indices[i]][lengthIndex];
			for (var j = 0; j < length; j++) {
				//Shift elements down by one (matlab 1-based index to 0-based index)
				connectedElement = md.mesh.vertexconnectivity[hit.indices[i]][j] - 1;
				indices = md.mesh.elements[connectedElement];
				connectedIndices.add(indices[0] - 1);
				connectedIndices.add(indices[1] - 1);
				connectedIndices.add(indices[2] - 1);
			}
		}

		//Apply modifications to included vertices in mesh using brush strength and falloff.
		var strength;
		for (var index of connectedIndices) {
			if (!baseIndices.has(index)) {
				strength = canvas.brush.strength * canvas.brush.falloff;
			}
			else {
				strength = canvas.brush.strength;
			}
			vertices[index*3+2] += strength * 100;
			md.geometry.surface[index] += strength;	
			md.geometry.thickness[index] += strength;
			coords[index*2+1] += strength;
			canvas.unitData[index] += strength;
		}
		
		//Update mesh on GPU
		bufferVertices.upload(canvas.gl.DYNAMIC_DRAW);
		bufferCoords.upload(canvas.gl.DYNAMIC_DRAW);
		canvas.unitNode.mesh.octree = new GL.Octree(canvas.unitNode.mesh);	
		
		//Update clouds if rendered
		//TODO: Steven, once you update the cloud generation in applyoptions.js, modify this code block to move the clouds as well. We'll want to move them individually later, but moving them all is ok for now.
		for (var i = 0; i < canvas.clouds.quantity; i++) {
			if (canvas.nodes['clouds' + i]) {
				var v1 = vec3.fromValues(vertices[hit.indices[0] * 3], vertices[hit.indices[0] * 3 + 1], vertices[hit.indices[0] * 3 + 2]);
				var v2 = vec3.fromValues(vertices[hit.indices[1] * 3], vertices[hit.indices[1] * 3 + 1], vertices[hit.indices[1] * 3 + 2]);
				var v3 = vec3.fromValues(vertices[hit.indices[2] * 3], vertices[hit.indices[2] * 3 + 1], vertices[hit.indices[2] * 3 + 2]);
				vec3.transformMat4(v1, v1, canvas.unitNode.modelMatrix);
				vec3.transformMat4(v2, v2, canvas.unitNode.modelMatrix);
				vec3.transformMat4(v3, v3, canvas.unitNode.modelMatrix);
				var x  = (v1[0] + v2[0] + v3[0]) / 3 + Math.floor((Math.random() * (1 + 10000 - (-10000)) + (-10000)));
				var y  = (v1[1] + v2[1] + v3[1]) / 3 + Math.floor((Math.random() * (1 + 10000 - (-10000)) + (-10000)));
				var z  = (v1[2] + v2[2] + v3[2]) / 3;
				canvas.nodes['clouds' + i].translation = [x, y + canvas.clouds.height, z];
				updateModelMatrix(canvas.nodes['clouds' + i]);
			}
		}
	}
} //}}}
function initializeMarker(canvas) { //{{{
	//Initialize data marker and tooltip display once per page load
	var marker = $('#' + canvas.dataMarkers.marker.id);
	var size = canvas.dataMarkers.size;
	if (!marker.hasClass('tooltipstered')) {
		marker.css({
			'position': 'absolute',
			'left': -size[0] + 'px',
			'top': -size[1] + '0px',
			'width': size[0] + 'px',
			'height': size[1] + 'px',
			'pointer-events': 'all',
			'cursor': 'pointer',
			'display': 'none'
		});
		marker.tooltipster({
			contentAsHTML: true,
			maxWidth: 320,
			maxHeight: 320,
			zIndex: 1000,
			trigger: 'custom',
			triggerOpen: {
				mouseenter: false,
				originClick: true,
				touchstart: false
			},
			triggerClose: {
				mouseleave: false,
				originClick: true,
				touchleave: false
			},
		});
		marker.on('click touch', function () {
			marker.fadeOut(175);
			marker.tooltipster('close');
		});
		canvas.dataMarkers.marker.selector = marker;
	}
	updateMarker(canvas, true);
} //}}}
function repositionMarker(canvas) { //{{{
	//Mover marker to point to mouse position, offset in y by 1 to enable immediate clicking.
	if (isEmptyOrUndefined(canvas.dataMarkers.marker.hit) || !canvas.camera.ready) { return; }
	var size = canvas.dataMarkers.size;
	var screenPoint = vec3.transformMat4(vec3.create(), canvas.dataMarkers.marker.hit.pos, canvas.camera.vpMatrix);
	//console.log(canvas, canvas.selector, $(canvas.id)
	var x = (screenPoint[0] + 1.0) * (canvas.width / 2) + canvas.selector.offset().left;
	var y = (-screenPoint[1] + 1.0) * (canvas.height / 2) + canvas.selector.offset().top;
	canvas.dataMarkers.marker.selector.css({
		'left': (Math.round(x) - size[0] / 2) + 'px', 
		'top': (Math.round(y) - size[1] + 1) + 'px'
	});
	
	if (canvas.dataMarkers.marker.selector.tooltipster('status').state != 'closed') { canvas.dataMarkers.marker.selector.tooltipster('reposition'); }
} //}}}
function updateMarker(canvas, reset) { //{{{
	//Retrieve data value fields and plots them on data marker popup if a hit has been registered.
	//TODO: Automatically pick up any field of size md.mesh.numberofelements
	//If no marker has been placed, no update is needed. If canvas is resimulating and unitNode has not been set yet, wait and try again.
	if (isEmptyOrUndefined(canvas.dataMarkers.marker.hit)) { return; }
	if (isEmptyOrUndefined(canvas.unitNode)) { setTimeout( function(){ updateMarker(canvas, reset); }, 750); return; }
	
	var hit = canvas.dataMarkers.marker.hit;
	
	var coords = canvas.unitNode.mesh.vertexBuffers.coords.data;
	var latitude = md.mesh.lat;
	var longitude = md.mesh.long;
	var thickness;
	var velocity;
	if (md.results[0]) {
		thickness = md.results[canvas.animation.frame].Thickness;
		velocity = md.results[canvas.animation.frame].Vel;
	}
	else {
		thickness = md.geometry.thickness;
		velocity = md.initialization.vel;
	}
	
	//Determine data values at hit position.
	var hitCoords = [coords[hit.indices[0]*2], coords[hit.indices[0]*2+1], coords[hit.indices[1]*2], coords[hit.indices[1]*2+1], coords[hit.indices[2]*2], coords[hit.indices[2]*2+1]];
	var hitLatitude = [latitude[hit.indices[0]], latitude[hit.indices[1]], latitude[hit.indices[2]]];
	var hitLongitude = [longitude[hit.indices[0]], longitude[hit.indices[1]], longitude[hit.indices[2]]];
	var hitThickness = [thickness[hit.indices[0]], thickness[hit.indices[1]], thickness[hit.indices[2]]];
	var hitVelocity = [velocity[hit.indices[0]], velocity[hit.indices[1]], velocity[hit.indices[2]]];
	var u = hitCoords[0] * hit.uvw[0] + hitCoords[2] * hit.uvw[1] + hitCoords[4] * hit.uvw[2];
	var v = hitCoords[1] * hit.uvw[0] + hitCoords[3] * hit.uvw[1] + hitCoords[5] * hit.uvw[2];
	var value = canvas.unitNode.caxis[0] * (1.0 - v) + canvas.unitNode.caxis[1] * v;
	var valueLatitude = Math.abs(hitLatitude[0] * hit.uvw[0] + hitLatitude[1] * hit.uvw[1] + hitLatitude[2] * hit.uvw[2]);
	var valueLongitude = Math.abs(hitLongitude[0] * hit.uvw[0] + hitLongitude[1] * hit.uvw[1] + hitLongitude[2] * hit.uvw[2]);
	var valueThickness = hitThickness[0] * hit.uvw[0] + hitThickness[1] * hit.uvw[1] + hitThickness[2] * hit.uvw[2];
	var valueVelocity = hitVelocity[0] * hit.uvw[0] + hitVelocity[1] * hit.uvw[1] + hitVelocity[2] * hit.uvw[2];	

	//Construct new argument array of the data display format for sprintf using first first argument as the formatSpecifier string and the rest as the additional arguments.
	var format = canvas.dataMarkers.format.slice();
	for (var i = 1; i < format.length; i++) {
		if (format[i].toLowerCase() == 'x') { format[i] = hit.modelPos[0]; }
		else if (format[i].toLowerCase() == 'y') { format[i] = hit.modelPos[1]; }
		else if (format[i].toLowerCase() == 'z') { format[i] = hit.modelPos[2]; }
		else if (format[i].toLowerCase() == 'lat') { format[i] = valueLatitude; }
		else if (format[i].toLowerCase() == 'long') { format[i] = valueLongitude; }
		else if (format[i].toLowerCase() == 'thickness') { format[i] = valueThickness; }
		else if (format[i].toLowerCase() == 'vel') { format[i] = valueVelocity; }
		else if (format[i].toLowerCase() == 'value') { format[i] = value; }
	}
	
	//Apply changes to tooltip
	$('#tooltip-content-data-marker-' + canvas.id).html(sprintf.apply(null, format));
	$('#tooltip-content-data-marker-' + canvas.id).css({'font': canvas.dataMarkers.font});				
	
	//If animated, setup animation loop to update plot as movie plays.
	if (canvas.dataMarkers.animated) {
		var isEmpty = (canvas.dataMarkers.values.length == 0);
		var lastUpdatedIndex = (canvas.dataMarkers.values.length-1);
		var newMovieFrame = (!isEmpty && canvas.dataMarkers.values[lastUpdatedIndex][0] != canvas.animation.frame);
		//If new data marker has been placed, reinitialize plot. If not, push new value into plot value array.
		if (reset) {
			canvas.dataMarkers.values = [];
			newMovieFrame = true;
			for (var currentFrame = 0; currentFrame < (canvas.unitNode.movieLength); currentFrame++) {
				coords = canvas.unitNode.texcoords[currentFrame];
				var hitCoords = [coords[hit.indices[0]*2], coords[hit.indices[0]*2+1], coords[hit.indices[1]*2], coords[hit.indices[1]*2+1], coords[hit.indices[2]*2], coords[hit.indices[2]*2+1]];
				var u = hitCoords[0] * hit.uvw[0] + hitCoords[2] * hit.uvw[1] + hitCoords[4] * hit.uvw[2];
				var v = hitCoords[1] * hit.uvw[0] + hitCoords[3] * hit.uvw[1] + hitCoords[5] * hit.uvw[2];
				var value = canvas.unitNode.caxis[0] * (1.0 - v) + canvas.unitNode.caxis[1] * v;
				canvas.dataMarkers.values.push([currentFrame, value]);
			}
		}
		else {
			if (isEmpty || newMovieFrame) {
				canvas.dataMarkers.values.push([canvas.animation.frame, value]);
			}
		}
		
		//Replot data marker popup using update data value array.
		if (isEmpty || newMovieFrame) {
			var dataLabels = {'latitude': valueLatitude, 'longitude': valueLongitude, 'thickness': valueThickness, 'velocity': valueVelocity, 'value': value};
			var dataDisplay = canvas.dataMarkers.values.slice(0, canvas.animation.frame+1);					
			plot(
				'id', '#sim-plot', 
				'type', 'bar', 
				'width', 400, 
				'height', 300, 
				'nticks', 25, 
				'xlabel', 'Time', 
				'ylabel', 'Value', 
				'title', 'Changes Over Time', 
				'datalabels', canvas.dataMarkers.labels,
				'labelvalues', dataLabels,
				'data', dataDisplay
			);
		}
	}
	repositionMarker(canvas);
	if (reset) {
		canvas.dataMarkers.marker.selector.fadeIn(175);
		canvas.dataMarkers.marker.selector.tooltipster('open');
	}
} //}}}
//}}}
//{{{ Drawing Functions
function updateCameraMatrix(canvas) { //{{{
    //Update view matrix and multiply with projection matrix to get the view-projection matrix.
	var vMatrix = mat4.create();
	var pMatrix = mat4.create();
	var translateMatrix = mat4.create();
	var rotationMatrix = mat4.create();
	var azimuthRotationMatrix = mat4.create();
	var elevationRotationMatrix = mat4.create();
	var aspectRatio = canvas.clientWidth / canvas.clientHeight;
	var cameraPosition = vec3.create();

	if (canvas.view.twod) { mat4.ortho(pMatrix, -aspectRatio*6.371e6/canvas.view.zoom, aspectRatio*6.371e6/canvas.view.zoom, -6.371e6/canvas.view.zoom, 6.371e6/canvas.view.zoom, canvas.camera.near, canvas.camera.far); }
	else { mat4.perspective(pMatrix, canvas.camera.fov * DEG2RAD, aspectRatio, canvas.camera.near, canvas.camera.far); }
	
	//Apply worldspace translation
	mat4.translate(vMatrix, translateMatrix, vec3.negate(vec3.create(), canvas.view.position));
	
	//Calculate rotation around camera focal point about worldspace origin
	if (canvas.view.twod) {
		mat4.rotate(azimuthRotationMatrix, azimuthRotationMatrix, DEG2RAD * 0, [0, 1, 0]);
		mat4.rotate(elevationRotationMatrix, elevationRotationMatrix, DEG2RAD * 90, [1, 0, 0]);
		mat4.multiply(rotationMatrix, elevationRotationMatrix, azimuthRotationMatrix);
	}
	else {
		mat4.rotate(azimuthRotationMatrix, azimuthRotationMatrix, DEG2RAD * canvas.view.rotation[0], [0, 1, 0]);
		mat4.rotate(elevationRotationMatrix, elevationRotationMatrix, DEG2RAD * canvas.view.rotation[1], [1, 0, 0]);
		mat4.multiply(rotationMatrix, elevationRotationMatrix, azimuthRotationMatrix);
	}

	//Apply rotation transform
	mat4.multiply(vMatrix, rotationMatrix, vMatrix);
	
	//Apply screenspace translation to emulate rotation around point
	mat4.identity(translateMatrix);
	mat4.translate(translateMatrix, translateMatrix, [0.0, 0.0, -6.371e6/canvas.view.zoom]);
	mat4.multiply(vMatrix, translateMatrix, vMatrix);
	
	//Apply projection matrix to get camera matrix
	mat4.copy(canvas.camera.vMatrix, vMatrix);
	mat4.multiply(canvas.camera.vpMatrix, pMatrix, vMatrix);
	
	//Calculate inverse view matrix fields for lighting and raycasts
	mat4.invert(canvas.camera.vInverseMatrix, canvas.camera.vMatrix);
	mat4.invert(canvas.camera.vpInverseMatrix, canvas.camera.vpMatrix);
	
	vec3.transformMat4(canvas.camera.position, cameraPosition, canvas.camera.vpInverseMatrix);
	canvas.camera.ready = true;
}//}}}
function drawSceneGraphNode(canvas, node) { //{{{
	if (!node.enabled) { return; }

	var gl = canvas.gl;
	gl.makeCurrent();
	
	var mvpMatrix = mat4.create();
	mat4.multiply(mvpMatrix, canvas.camera.vpMatrix, node.modelMatrix);
	
	var mvMatrix = mat4.create();
	mat4.multiply(mvMatrix, canvas.camera.vMatrix, node.modelMatrix);
	
	var normalMatrix = mat4.create();
	mat4.invert(normalMatrix, mvMatrix);
	mat4.transpose(normalMatrix, normalMatrix);
	
	if (node.texture) { node.texture.bind(0); }
	if (node.disableDepthTest) { gl.disable(gl.DEPTH_TEST); }
	if (node.enableCullFace) { gl.enable(gl.CULL_FACE); }

	gl.cullFace(node.cullFace);
	gl.lineWidth(node.lineWidth);
	gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

	//Setup for light that originates from camera
	var origin = vec3.fromValues(0, 0, 0);
	var lightOrigin = vec3.fromValues(0, 0, 0);
	var cameraPositionRelative = vec3.create();
	vec3.transformMat4(origin, origin, canvas.camera.vInverseMatrix);
	vec3.normalize(lightOrigin, lightOrigin);
	vec3.sub(cameraPositionRelative, origin, node.translation);
	cameraHeight = vec3.length(cameraPositionRelative);
	
	var atm = { 							//Default Values
				wavelength_r: 0.65, 		//0.65		Red wavelength (micrometers)
				wavelength_g: 0.57,			//0.57		Green wavelength (micrometers)
				wavelength_b: 0.475,		//0.475		Green wavelength (micrometers)
				eSun: 100.0,				//20.0		Sun intensity	
				kRayleigh: 0.0025,			//0.0025	Rayleigh scattering amount
				kMie: 0.000, 				//0.01		Mie scattering amount
				g: -0.99,					//-0.99		Mie phase asymmetry/direction factor
				hdr_exposure: 0.8,			//0.8		High Dynamic Range Exposure
				scale: 1.25, 				//1.025		Scale of atmosphere. WARNING: Change atmosphereScale in applyoptions.js, and scaling constants.
				scaleDepth: 0.25, 			//0.25		Percentage altitude at which the atmosphere's average density is found
				a: -0.00287,				//-0.00287	Scaling constant a
				b: 0.459,					//0.459		Scaling constant b
				c: 3.83,					//3.83		Scaling constant c
				d: -6.80,					//-6.80		Scaling constant d
				e: 3.6,						//5.25		Scaling constant e. Lower when increasing atmosphere scale.
				attenuation: 0.5			//0.5		Strength of atmospheric scattering on ground shading.
	};
			
	var inv_wavelength4 = [1.0 / Math.pow(atm.wavelength_r, 4), 1.0 / Math.pow(atm.wavelength_g, 4), 1.0 / Math.pow(atm.wavelength_b, 4)];
	var innerRadius = 6.371e6;
	var outerRadius = innerRadius*atm.scale;
	var scale = 1.0 / (outerRadius - innerRadius);
	var scaleDepth = atm.scaleDepth;
	
	node.shader.uniforms({
		m4MVP: mvpMatrix,
		m4Normal: normalMatrix,
		m4Model: node.modelMatrix,
		//u_lightPosition: [-lightOrigin[0], -lightOrigin[1], -lightOrigin[2]],
		u_lightPosition: [1.0, 1.0, 1.0],
		u_diffuseColor: [1.0, 0.9, 0.9],
		u_texture: 0,
		u_alpha: node.alpha,
		u_maskEnabled: node.maskEnabled,
		u_maskHeight: node.maskHeight,
		u_maskColor: node.maskColor,
		v3CameraPosition: origin,
		v3Translate: node.translation,
		v3LightPos: lightOrigin,
		v3InvWavelength: inv_wavelength4,
		fOuterRadius: outerRadius,
		fOuterRadius2: outerRadius * outerRadius,
		fInnerRadius: innerRadius,
		fInnerRadius2: innerRadius * innerRadius,
		fKrESun: atm.kRayleigh * atm.eSun, 
		fKmESun: atm.kMie * atm.eSun, 
		fKr4PI: atm.kRayleigh * 4 * Math.PI, 
		fKm4PI: atm.kMie * 4 * Math.PI,
		fScale: scale, 
		fScaleDepth: scaleDepth,
		fScaleOverScaleDepth: scale/scaleDepth, 
		v3LightPosFrag: lightOrigin,
		fHdrExposure: atm.hdr_exposure,	
		g: atm.g,			
		g2: atm.g * atm.g,
		a: atm.a,
		b: atm.b,
		c: atm.c,
		d: atm.d,		
		e: atm.e,
		attenuation: atm.attenuation
	}).draw(node.mesh, node.drawMode, 'triangles');
	
	gl.enable(gl.DEPTH_TEST);
	gl.disable(gl.CULL_FACE);
} //}}}
function draw(canvas) { //{{{
	//Ensure all nodes are ready to render
	//TODO: Come up with better way to check if shaders are ready, or move outside of main draw function
	var nodes = canvas.nodes;
	if (!canvas.draw.ready) {
		if (nodes.length !== 0) {
			canvas.draw.ready = true;
			for (var node in nodes) {
				if (nodes[node].shader.ready == false) {
					canvas.draw.ready = false;
					break;
				}
			}
			
		}
	}
	
	//Begin rendering nodes
	if (canvas.draw.ready) {
		if (canvas.textcanvas) { canvas.textcanvas.draw(canvas); }
		if (canvas.overlaycanvas) { canvas.overlaycanvas.draw(canvas); }
	
		var rect = canvas.getBoundingClientRect();
		canvas.width  = rect.width;
		canvas.height = rect.height;
		
		var gl = canvas.gl;
		gl.makeCurrent(); //litegl function to handle switching between multiple canvases
		gl.viewport(0, 0, canvas.width, canvas.height);
		gl.clearColor(canvas.backgroundcolor[0], canvas.backgroundcolor[1], canvas.backgroundcolor[2], canvas.backgroundcolor[3]);
		gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
		
		updateCameraMatrix(canvas);
		
		var drawPassNumber = 3;
		for (var i = drawPassNumber - 1; i >= 0; i--) {
			for (var node in nodes) {
				if (nodes[node].drawOrder == i) { drawSceneGraphNode(canvas, nodes[node]); }
			}
		}
	}
	
	//Regardless of ready state, schedule next frame to check for ready state and render
	canvas.draw.handler = window.requestAnimationFrame(function(time) { draw(canvas); });
} //}}}
//}}}
