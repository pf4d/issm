//SETTINGS class definition
//
//   Usage:
//      settings=new settings();

function settings (){
	//methods
	this.setdefaultparameters = function(){// {{{
		//are we short in memory ? (0 faster but requires more memory)
		this.lowmem=0;

		//i/o:
		this.io_gather=1;

		//results frequency by default every step
		this.output_frequency=1;

		//checkpoints frequency, by default never: 
		this.recording_frequency=0;

		//this option can be activated to load automatically the results
		//onto the model after a parallel run by waiting for the lock file
		//N minutes that is generated once the solution has converged
		//0 to deactivate
		this.waitonlock=Infinity;

		//upload options: 
		upload_port         = 0;
		
		//throw an error if solver residue exceeds this value
		self.solver_residue_threshold=1e-6;

	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   settings class echo:'));
		
		fielddisplay(this,'results_on_nodes','results are output for all the nodes of each element');
		fielddisplay(this,'io_gather','I/O gathering strategy for result outputs (default 1)');
		fielddisplay(this,'lowmem','is the memory limited ? (0 or 1)');
		fielddisplay(this,'output_frequency','frequency at which results are saved in all solutions with multiple time_steps');
		fielddisplay(this,'recording_frequency','frequency at which the runs are being recorded, allowing for a restart');
		fielddisplay(this,'waitonlock','maximum number of minutes to wait for batch results (NaN to deactivate)');
		fielddisplay(this,'upload_server','server hostname where model should be uploaded');
		fielddisplay(this,'upload_path','path on server where model should be uploaded');
		fielddisplay(this,'upload_login','server login');
		fielddisplay(this,'upload_port','port login (default is 0)');
		fielddisplay(this,'upload_filename','unique id generated when uploading the file to server');
		fielddisplay(this,'solver_residue_threshold','throw an error if solver residue exceeds this value');


	}// }}}
	this.classname= function(){// {{{
		return "settings";

	}// }}}
		this.checkconsistency = function(md,solution,analyses) { // {{{

			checkfield(md,'fieldname','settings.results_on_nodes','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','settings.io_gather','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','settings.lowmem','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','settings.output_frequency','numel',[1],'>=',1);
			checkfield(md,'fieldname','settings.recording_frequency','numel',[1],'>=',0);
			checkfield(md,'fieldname','settings.waitonlock','numel',[1]);
			checkfield(md,'fieldname','settings.solver_residue_threshold','numel',[1],'>',0);
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'object',this,'fieldname','results_on_nodes','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','io_gather','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','lowmem','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','output_frequency','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','recording_frequency','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','solver_residue_threshold','format','Double');
			if (this.waitonlock>0) WriteData(fid,prefix,'name','md.settings.waitonlock','data',true,'format','Boolean');
			else WriteData(fid,prefix,'name','md.settings.waitonlock','data',false,'format','Boolean');
		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{
	this.results_on_nodes    = 0;
	this.io_gather           = 0;
	this.lowmem              = 0;
	this.output_frequency    = 0;
	this.recording_frequency   = 0;
	this.waitonlock          = 0;
	this.upload_server       = '';
	this.upload_path         = '';
	this.upload_login        = '';
	this.upload_port         = 0;
	this.upload_filename     = '';
	this.solver_residue_threshold = 0;
	this.setdefaultparameters();
	//}}}
}
