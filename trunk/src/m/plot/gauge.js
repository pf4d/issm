function gaugeInit(){ //{{{
	/*
		References:	https://bernii.github.io/gauge.js/
	*/
	// Convert arguments to options.
	var args 		= Array.prototype.slice.call(arguments);
	var options 	= new pairoptions(args);

	// Recover option values.
	var canvasId	= options.getfieldvalue('canvasId', 'sim-gauge-canvas');
	var textFieldId	= options.getfieldvalue('textFieldId', 'gauge-text');
	var	value 		= options.getfieldvalue('value', 0);
	var max 		= options.getfieldvalue('max', 22);
	var colors 		= options.getfieldvalue('colors', [[0.0, "#FFFF00"],
													   [0.50, "#FF8000"],
													   [1.0, "#FF0000"]]);
	
	// Construct associative array of options.
	var opts = {
		lines: 				12, 	// The number of lines to draw
		angle: 				0,		// The length of each line
		lineWidth: 			0.39,	// The line thickness
		
		pointer: {
			length: 			0, 			// The radius of the inner circle
			strokeWidth: 		0, 			// The rotation offset
			color: 				'#000000' 	// Fill color
		},
		
		limitMax: 			'false',	// If true, the pointer will not go past the end of the gauge
		percentColors: 		colors,
		strokeColor: 		'#E0E0E0',	// to see which ones work best for you
		generateGradient: 	true
	};
	
	// Create new gauge.
	var newGauge = new Gauge(document.getElementById(canvasId)).setOptions(opts);
	
	// Bind text field to gauge.
	newGauge.setTextField(document.getElementById(textFieldId));
	
	newGauge.textField.render = function(newGauge){
    	return this.el.innerHTML = newGauge.displayedValue.toFixed(0);
    };
    
	newGauge.maxValue 			= max;
	newGauge.animationSpeed 	= 1;
	newGauge.set(1);		// Somehow this prevents pointer from rendering on initialization.
	newGauge.set(value);	// set actual value
	
	return newGauge;
} //}}}
