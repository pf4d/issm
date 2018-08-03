function plotmodel(md){ //{{{

	//Convert arguments to array: 
	var args = Array.prototype.slice.call(arguments);

	//First process options
	var  options = new plotoptions(args.slice(1,args.length));

	
	//get number of subplots
	subplotwidth=Math.ceil(Math.sqrt(options.numberofplots)); 
	
	//Get figure number and number of plots
	numberofplots=options.numberofplots;

	//if nlines and ncols specified, then bypass.
	var nlines,ncols;
	if (options.list[0].exist('nlines')){
		nlines=options.list[0].getfieldvalue('nlines');
	}
	else {
		nlines=Math.ceil(numberofplots/subplotwidth);
	}
	if (options.list[0].exist('ncols')){
		ncols=options.list[0].getfieldvalue('ncols');
	}
	else {
		ncols=subplotwidth;
	}
	
	//check that nlines and ncols were given at the same time!
	if ((options.list[0].exist('ncols') & !options.list[0].exist('nlines')) | (options.list[0].exist('nlines') & !options.list[0].exist('ncols'))) throw Error('plotmodel error message: nlines and ncols  need to be specified together, or not at all');

	//go through subplots
	if (numberofplots){
		//Reinitialize all canvases
		for (var i=0;i<numberofplots;i++){
			document.getElementById(options.list[i].getfieldvalue('canvasid')).initialized = false;
		}
		//Go through all data plottable and close window if an error occurs
		for (var i=0;i<numberofplots;i++){
			plot_manager(options.list[i].getfieldvalue('model',md),options.list[i],subplotwidth,nlines,ncols,i);

			//List all unused options
			options.list[i].displayunused();
		}
	}
} //}}}
