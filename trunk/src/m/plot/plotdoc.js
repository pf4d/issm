function plotdoc() { //{{{
	//PLOTDOC - plot documentation
	//
	//   Usage:
	//      plotdoc()
	
	//TODO: Standardize image to overlay_image, heightscale to scaling, colorbarfontsize/color, clarify innermask/outermask, edgecolor implementation, check colormap, 

	console.log(' WARNING: starred methods (*) are experimental and not guarenteed to be stable');
	console.log('   Plot usage: plotmodel(model,varargin)');
	console.log('   Options: ');
	console.log('       "canvasid": canvas id');
	console.log('       "data" : what we want to plot');
	console.log('                Available values for "data" are: ');
	console.log('                  - any field of the model structure. ex: plot(md,"data","vel"), or plot(md,"data",md.initialization.vel)');
	console.log('                  - "mesh": draw mesh using trisurf');
	console.log('                  - "quiver": quiver plot');
	console.log('       "2d": renders orthographic camera with view set to [0, 90] (default "off", ex: "on", "off")');
	console.log('       "backgroundcolor": plot background color. (default "lightcyan", ex: "green","blue")');
	console.log('       "brush": specify brush options (default {"strength":0.075,"falloff":0.5})');
	console.log('       "caxis": modify  colorbar range. (array of type [a, b] where b>=a)');
	console.log('       "colorbar": add colorbar (default "off", ex: "on", "off")');
	console.log('       "colorbarid": colorbar canvas id (string)');
	console.log('       "colorbartitle": colorbar title (string)');	
	console.log('       "colorbarnticks": number of colorbar ticks (default 6, ex: 2, 10)');	
	console.log('       "colorbarprecision": colorbar label digit precision (default 3, ex: 1, 4)');	
	console.log('       "colorbarinnerlabels": choose whether labels are inside colorbar (default "off", ex: "on", "off")');
	console.log('       "colorbarfontsize": specify colorbar font size (default 1, ex: 14, 22)');
	console.log('       "colorbarfontcolor": specify colorbar font color (default "black", ex: "green","blue")');
	console.log('       "colorbarwidth": multiplier (default 1) to the default width colorbar');
	console.log('       "colorbarheight": multiplier (default 1) to the default height colorbar');
	console.log('       "colormap": same as standard matlab option (default "jet", ex: "hsv","cool","spring","gray","Ala","Rignot",...)');
	console.log('       "controlsensitivity": sensitivty of view/zoom changes as a percentage of default (default 1, ex: 0.5, 2.75)');
	console.log('       "datamarkers": toggle data marker displays (default "on", ex: "on", "off")');
	console.log('       "datamarkers_image": toggle data marker displays (default "on", ex: "on", "off")');
	console.log('       "datamarkerssize": specifiy the width and height of the data markers (default [32,32], ex: [24,32], [105,10])');
	console.log('       "datamarkersoptions": specifiy options for data markers (default {"enabled":"on","image":canvas.rootPath+"textures/data_marker.svg","size":[32,32],"format":["X: %.2e<br>Y: %.2e<br>Z: %.2e]<br>Value: %0.1f","x","y","z","value"]}');
	console.log('       "displayview": print view value to console (default "off", ex: "on", "off")');
	console.log('       "displayzoom": print zoom value to console (default "off", ex: "on", "off")');
	console.log('       "edgecolor": same as standard matlab option EdgeColor (default "black", ex: color name: "blue" or RGB array: [0.5, 0.2, 0.8])');
	console.log('       "heightscale": scaling factor to accentuate height. (default 1, ex: 0.5, 100)');
	console.log('       "linewidth*": line width for mesh, quiver, and contour plots, currently limited by WebGL to 1. (default 1, ex: 2, 5)');
	console.log('       "log": value of log (default 10, ex: 2, Math.E)');
	console.log('       "mask": list of flags of size numberofnodes or numberofelements. Only "true" values are plotted ');
	console.log('       "movieoptions": specify movie options (default {"fps":4,"loop":true})');
	console.log('       "innermask*": Special mask that colors all parts of a data mesh below a height a certain color. provide innermaskheight and innermaskcolor options also (default "off", ex: "on", "off")');
	console.log('       "outermask*": Special mask that colors all parts of a overlay mesh below a height a certain color. provide outermaskheight and outermaskcolor options also (default "off", ex: "on", "off")');
	console.log('       "overlay": overlay a radar amplitude image behind (default "off", ex: "on", "off")');
	console.log('       "overlay_image": path to overlay image (default "", ex: "./images/radar.png")');
	console.log('       "quiver": add quiver plot overlay for velocities. (default "off", ex: "on", "off")');
	console.log('       "scaling": scaling factor used by quiver plots. Default is 0.4');
	console.log('       "alpha": transparency coefficient 0.0 to 1.0, the lower, the more transparent. (default 1.0, ex: 0.5, 0.25)');
	console.log('       "azlim": azimuth view limits (ex: [0, 180])');
	console.log('       "ellim": elevation view limits (ex: [-90, 90])');
	console.log('       "origin": initial camera offset from model center (default [0,0,0.0], ex: [-2, 1.5, 0.01])');
	console.log('       "render": toggle sky, ground, and space rendering. (default [], ex: ["sky", "space"], ["ground"])');
	console.log('       "viewPanning": enable view origin panning with two-finger touch or shift+mouse drag. (default "off", ex: "on", "off")');
	console.log('       "view": initial azimuth and elevation angles for camera (default [0,90], ex: [90, 180]');
	console.log('       "xlim": x coordinates to fit inside camera (ex: [0, 500])');
	console.log('       "ylim": y coordinates to fit inside camera (ex: [0, 500])');
	console.log('       "zlim": z coordinates to fit inside camera (ex: [0, 500])');
	console.log('       "zoomlim": zoom view limits (ex: [0.05, 10])');
	console.log('       "zoom": initial camera zoom as a percentage of default (default 1, ex: 1.5, 0.01)');
	console.log('       "cloud*": plot a cloud of points, given a flat array of 3d coordinates (ex: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0])');
	console.log('       "expdisp*": plot exp file on top of a data plot. provide exp file as an argument (use a cell of strings if more than one)');
	console.log('       "textlabels*": plot text labels rendered in 3d space, using an array of text/coordinate pairs (ex: [{"pos":[0.0,0.0,0.0],"text":"origin"}])');
	
	console.log('  ');
	console.log('   Examples:');
	console.log('       plotmodel(md,"data","vel","data","mesh","view#2",3,"colorbar#all","on")');
	console.log('       plotmodel(md,"data",md.geomtery.surface)');
} //}}}
