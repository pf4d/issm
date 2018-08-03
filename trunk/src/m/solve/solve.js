function solve(md,solutionstring){
//SOLVE - apply solution sequence for this model
//
//   Usage:
//      solve(md,solutionstring,varargin)
//      where varargin is a lit of paired arguments of string OR enums
//
//   solution types available comprise:
//		 - 'Stressbalance'      or 'sb'
//		 - 'Masstransport'      or 'mt'
//		 - 'Thermal'            or 'th'
//		 - 'Steadystate'        or 'ss'
//		 - 'Transient'          or 'tr'
//		 - 'Balancethickness'   or 'mc'
//      - 'Balancevelocity'   or 'bv'
//		 - 'BedSlope'           or 'bsl'
//		 - 'SurfaceSlope'       or 'ssl'
//		 - 'Hydrology'          or 'hy'
//      - 'DamageEvolution'   or 'da'
//		 - 'Gia'                or 'gia'
//		 - 'Sealevelrise'       or 'slr'
//
//  extra options:
//      - loadonly    : does not solve. only load results
//      - runtimename : true or false (default is true), makes name unique
//      - checkconsistency : 'yes' or 'no' (default is 'yes'), ensures checks on consistency of model
//      - restart: 'directory name (relative to the execution directory) where the restart file is located.
//      - callback: callback function to be called upon receiving the results from the server, or local computations. 
//
//   Examples:
//      md=solve(md,'Stressbalance');
//      md=solve(md,'sb');

	if(typeof solutionstring !== 'string') {
		throw Error(sprintf("%s\n", "ISSM's solve function only accepts strings for solution sequences. Type help solve to get a list of supported solutions."));
	}

	//recover and process solve options
	if((solutionstring.toLowerCase() === 'sb') || (solutionstring.toLowerCase() === 'stressbalance')){
		solutionstring = 'StressbalanceSolution';
	}else if((solutionstring.toLowerCase() === 'mt') || (solutionstring.toLowerCase() === 'masstransport')){
		solutionstring = 'MasstransportSolution';	
	}else if((solutionstring.toLowerCase() === 'th') || (solutionstring.toLowerCase() === 'thermal')){
		solutionstring = 'ThermalSolution';
	}else if((solutionstring.toLowerCase() === 'st') || (solutionstring.toLowerCase() === 'steadystate')){
		solutionstring = 'SteadystateSolution';
	}else if((solutionstring.toLowerCase() === 'tr') || (solutionstring.toLowerCase() === 'transient')){
		solutionstring = 'TransientSolution';
	}else if((solutionstring.toLowerCase() === 'mc') || (solutionstring.toLowerCase() === 'balancethickness')){
		solutionstring = 'BalancethicknessSolution';
	}else if((solutionstring.toLowerCase() === 'bv') || (solutionstring.toLowerCase() === 'balancevelocity')){
		solutionstring = 'BalancevelocitySolution';
	}else if((solutionstring.toLowerCase() === 'bsl') || (solutionstring.toLowerCase() === 'bedslope')){
		solutionstring = 'BedSlopeSolution';
	}else if((solutionstring.toLowerCase() === 'ssl') || (solutionstring.toLowerCase() === 'surfaceslope')){
		solutionstring = 'SurfaceSlopeSolution';
	}else if((solutionstring.toLowerCase() === 'hy') || (solutionstring.toLowerCase() === 'hydrology')){
		solutionstring = 'HydrologySolution';
	}else if((solutionstring.toLowerCase() === 'da') || (solutionstring.toLowerCase() === 'damageevolution')){
		solutionstring = 'DamageEvolutionSolution';
	}else if((solutionstring.toLowerCase() === 'gia') || (solutionstring.toLowerCase() === 'gia')){
		solutionstring = 'GiaSolution';
	}else if((solutionstring.toLowerCase() === 'slr') || (solutionstring.toLowerCase() === 'sealevelrise')){
		solutionstring = 'SealevelriseSolution';
	}else{
		throw Error(sprintf("%s%s%s\n",'solutionstring ',solutionstring,' not supported!'));
	}
	
	//Process options
	var args = Array.prototype.slice.call(arguments);
	var options = new pairoptions(args.slice(2,args.length));
	options.addfield('solutionstring',solutionstring);

	//recover some fields
	md.priv.solution=solutionstring;
	cluster=md.cluster;

	//check model consistency
	if (options.getfieldvalue('checkconsistency','yes') == 'yes'){
		if (md.verbose.solution){
			console.log('checking model consistency');
		}
		ismodelselfconsistent(md);
	}

	//If we are restarting, actually use the provided runtime name:
	restart=options.getfieldvalue('restart','');

	//First, build a runtime name that is unique
	if (restart==1 ){
		//Leave the runtimename as is
		}
	else{
		if (!(restart == '')){
			md.priv.runtimename=restart;
		}
		else if (options.getfieldvalue('runtimename',true)){
			c=new Date().getTime();
			md.priv.runtimename=sprintf('%s-%g',md.miscellaneous.name,c);
		}
		else{
			md.priv.runtimename=md.miscellaneous.name;
		}
	}

	//if running qmu analysis, some preprocessing of dakota files using models
	//fields needs to be carried out. 
	if (md.qmu.isdakota){
		throw Error("solve error message: qmu runs not supported yet!");
		//md.preqmu(options);
	}


	//Do we load results only?
	if (options.getfieldvalue('loadonly',false)){
		loadresultsfromcluster(md);
		return;
	}

	//Marshall into a binary array (fid) all the fields of model.
	var fid = marshall(md);                                          // bin file
	
	//deal with toolkits options: 
	toolkitsstring= md.toolkits.ToolkitsFile(md.miscellaneous.name + '.toolkits'); // toolkits file

	//callback function: 
	function callbackfunction(){solving=false;}; //default, do nothing if no callback function requested.
	if (options.getfieldvalue('callbackfunction',false)){
		callbackfunction=options.getfieldvalue('callbackfunction');
	}
	
	//callback error function: 
	function callbackerrorfunction(){solving=false;}; //default, do nothing if no callback function requested.
	if (options.getfieldvalue('callbackerrorfunction',false)){
		callbackerrorfunction=options.getfieldvalue('callbackerrorfunction');
	}
	
	//callback id: 
	var callbackid = '.run-button'; //default, update .run-button elements with progress updates.
	if (options.getfieldvalue('callbackid',false)){
		callbackid=options.getfieldvalue('callbackid');
	}

	if (cluster.classname() == 'local'){  //{{{

		/*We are running locally on the machine, using the issm module:*/
		console.log('running issm locally');
		
		//Call issm:
		var outputs = issm(fid, toolkitsstring, solutionstring, md.miscellaneous.name); 
		
		//Recover output arguments: 
		var outputbuffer = outputs[0]; var outputbuffersize = outputs[1];
			
		//Load results: 
		md = loadresultsfrombuffer(md,outputbuffer,outputbuffersize); 
		
		//Call back? 
		callbackfunction();

		return md;

	} //}}}
	else { //{{{

		/*We are running somewhere else on a computational server. Send the buffer to that server and retrieve output: */
		cluster.UploadAndRun(md,callbackfunction,callbackerrorfunction,callbackid,fid,toolkitsstring,solutionstring,md.miscellaneous.name,md.priv.runtimename);

		return md;

	} //}}}
}
