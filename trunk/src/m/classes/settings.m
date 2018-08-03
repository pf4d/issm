%SETTINGS class definition
%
%   Usage:
%      settings=settings();

classdef settings
	properties (SetAccess=public) 
		results_on_nodes    = 0;
		io_gather           = 0;
		lowmem              = 0;
		output_frequency    = 0;
		recording_frequency   = 0;
		waitonlock          = 0;
		upload_server       = '';
		upload_path         = '';
		upload_login        = '';
		upload_port         = 0;
		upload_filename     = '';
		solver_residue_threshold = 0;
	end
	methods
		function self = settings(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%are we short in memory ? (0 faster but requires more memory)
			self.lowmem=0;

			%i/o:
			self.io_gather=1;

			%results frequency by default every step
			self.output_frequency=1;

			%checkpoints frequency, by default never: 
			self.recording_frequency=0;

			%this option can be activated to load automatically the results
			%onto the model after a parallel run by waiting for the lock file
			%N minutes that is generated once the solution has converged
			%0 to deactivate
			self.waitonlock=Inf;

			%upload options: 
			self.upload_port         = 0;

			%throw an error if solver residue exceeds this value
			self.solver_residue_threshold = 1e-6;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','settings.results_on_nodes','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','settings.io_gather','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','settings.lowmem','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','settings.output_frequency','numel',[1],'>=',1);
			md = checkfield(md,'fieldname','settings.recording_frequency','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','settings.waitonlock','numel',[1]);
			md = checkfield(md,'fieldname','settings.solver_residue_threshold','numel',[1],'>',0);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   general settings parameters:'));

			fielddisplay(self,'results_on_nodes','results are output for all the nodes of each element');
			fielddisplay(self,'io_gather','I/O gathering strategy for result outputs (default 1)');
			fielddisplay(self,'lowmem','is the memory limited ? (0 or 1)');
			fielddisplay(self,'output_frequency','frequency at which results are saved in all solutions with multiple time_steps');
			fielddisplay(self,'recording_frequency','frequency at which the runs are being recorded, allowing for a restart');
			fielddisplay(self,'waitonlock','maximum number of minutes to wait for batch results (NaN to deactivate)');
			fielddisplay(self,'upload_server','server hostname where model should be uploaded');
			fielddisplay(self,'upload_path','path on server where model should be uploaded');
			fielddisplay(self,'upload_login','server login');
			fielddisplay(self,'upload_port','port login (default is 0)');
			fielddisplay(self,'upload_filename','unique id generated when uploading the file to server');
			fielddisplay(self,'solver_residue_threshold','throw an error if solver residue exceeds this value (NaN to deactivate)');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','results_on_nodes','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','io_gather','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','lowmem','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','output_frequency','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','recording_frequency','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','waitonlock','data',self.waitonlock>0,'format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','solver_residue_threshold','format','Double');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.settings.results_on_nodes'],self.results_on_nodes);
			writejsdouble(fid,[modelname '.settings.io_gather'],self.io_gather);
			writejsdouble(fid,[modelname '.settings.lowmem'],self.lowmem);
			writejsdouble(fid,[modelname '.settings.output_frequency'],self.output_frequency);
			writejsdouble(fid,[modelname '.settings.recording_frequency'],self.recording_frequency);
			writejsdouble(fid,[modelname '.settings.waitonlock'],self.waitonlock);
			writejsstring(fid,[modelname '.settings.upload_server'],self.upload_server);
			writejsstring(fid,[modelname '.settings.upload_path'],self.upload_path);
			writejsstring(fid,[modelname '.settings.upload_login'],self.upload_login);
			writejsdouble(fid,[modelname '.settings.upload_port'],self.upload_port);
			writejsstring(fid,[modelname '.settings.upload_filename'],self.upload_filename);
			writejsstring(fid,[modelname '.settings.solver_residue_threshold'],self.solver_residue_threshold);
		end % }}}
	end
end
