%QMU class definition
%
%   Usage:
%      qmu=qmu();

classdef qmu
	properties (SetAccess=public) 
		isdakota                    = 0;
		variables                   = struct();
		responses                   = struct();
		method                      = struct();
		params                      = struct();
		results                     = struct();
		partition                   = NaN;
		numberofpartitions          = 0;
		numberofresponses           = 0;
		variabledescriptors         = {};
		responsedescriptors         = {};
		mass_flux_profile_directory = NaN;
		mass_flux_profiles          = NaN;
		mass_flux_segments          = {};
		adjacency                   = NaN;
		vertex_weight               = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.partition=project3d(md,'vector',self.partition','type','node');
		end % }}}
		function self = qmu(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~md.qmu.isdakota, return; end

			version=IssmConfig('_DAKOTA_VERSION_'); version=str2num(version(1:3));

			if version < 6,
				if md.qmu.params.evaluation_concurrency~=1,
					md = checkmessage(md,['concurrency should be set to 1 when running dakota in library mode']);
				end
			else
				if ~strcmpi(self.params.evaluation_scheduling,'master'),
					md = checkmessage(md,['evaluation_scheduling in qmu.params should be set to ''master''']);
				end
				if md.cluster.np<=1,
					md = checkmessage(md,['in parallel library mode, Dakota needs to run on at least 2 cpus, 1 cpu for the master, 1 cpu for the slave. Modify md.cluser.np accordingly.']);
				end
					
				if self.params.processors_per_evaluation<1,
					md = checkmessage(md,['in parallel library mode, Dakota needs to run at least one slave on one cpu (md.qmu.params.processors_per_evaluation >=1)!']);
				end
				if mod(md.cluster.np-1,self.params.processors_per_evaluation), 
					md = checkmessage(md,['in parallel library mode, the requirement is for md.cluster.np = md.qmu.params.processors_per_evaluation * number_of_slaves, where number_of_slaves will automatically be determined by Dakota. Modify md.cluster.np accordingly']);
				end
			end
			if ~isempty(md.qmu.partition),
				if numel(md.qmu.partition)~=md.mesh.numberofvertices,
					md = checkmessage(md,['user supplied partition for qmu analysis should have size md.mesh.numberofvertices x 1 ']);
				end
				if min(md.qmu.partition)~=0,
					md = checkmessage(md,['partition vector not indexed from 0 on']);
				end
				if max(md.qmu.partition)>=md.qmu.numberofpartitions,
					md = checkmessage(md,['for qmu analysis, partitioning vector cannot go over npart, number of partition areas']);
				end
			end

			if ~strcmpi(md.cluster.name,'none'),
				if md.settings.waitonlock==0,
					md = checkmessage(md,['waitonlock should be activated when running qmu in parallel mode!']);
				end
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   qmu parameters:'));

			fielddisplay(self,'isdakota','is qmu analysis activated?');
			for i=1:numel(self.variables)
				disp(sprintf('         variables%s:  (arrays of each variable class)',...
					string_dim(self.variables,i)));
				fnames=fieldnames(self.variables(i));
				maxlen=0;
				for j=1:numel(fnames)
					maxlen=max(maxlen,length(fnames{j}));
				end

				for j=1:numel(fnames)
					disp(sprintf(['            %-' num2str(maxlen+1) 's:    [%ix%i]    ''%s'''],...
						fnames{j},size(self.variables.(fnames{j})),class(self.variables.(fnames{j}))));
				end
			end
			for i=1:numel(self.responses)
				disp(sprintf('         responses%s:  (arrays of each response class)',...
					string_dim(self.responses,i)));
				fnames=fieldnames(self.responses(i));
				maxlen=0;
				for j=1:numel(fnames)
					maxlen=max(maxlen,length(fnames{j}));
				end

				for j=1:numel(fnames)
					disp(sprintf(['            %-' num2str(maxlen+1) 's:    [%ix%i]    ''%s'''],...
						fnames{j},size(self.responses.(fnames{j})),class(self.responses.(fnames{j}))));
				end
			end
			fielddisplay(self,'numberofresponses','number of responses') 
			for i=1:numel(self.method);
				if strcmp(class(self.method(i)),'dakota_method')
					disp(sprintf('            method%s :    ''%s''',...
						string_dim(self.method,i),self.method(i).method));
				end
			end
			for i=1:numel(self.params)
				disp(sprintf('         params%s:  (array of method-independent parameters)',...
					string_dim(self.params,i)));
				fnames=fieldnames(self.params(i));
				maxlen=0;
				for j=1:numel(fnames)
					maxlen=max(maxlen,length(fnames{j}));
				end

				for j=1:numel(fnames)
					disp(sprintf(['            %-' num2str(maxlen+1) 's: %s'],...
						fnames{j},any2str(self.params(i).(fnames{j}))));
				end
			end
			for i=1:numel(self.results)
				disp(sprintf('         results%s:  (information from dakota files)',...
					string_dim(self.results,i)));
				fnames=fieldnames(self.results(i));
				maxlen=0;
				for j=1:numel(fnames)
					maxlen=max(maxlen,length(fnames{j}));
				end

				for j=1:numel(fnames)
					disp(sprintf(['            %-' num2str(maxlen+1) 's:    [%ix%i]    ''%s'''],...
						fnames{j},size(self.results.(fnames{j})),class(self.results.(fnames{j}))));
				end
			end
			fielddisplay(self,'partition','user provided mesh partitioning, defaults to metis if not specified') 
			fielddisplay(self,'numberofpartitions','number of partitions for semi-discrete qmu') 
			fielddisplay(self,'variabledescriptors','');
			fielddisplay(self,'responsedescriptors','');
			fielddisplay(self,'method','array of dakota_method class');
			fielddisplay(self,'mass_flux_profile_directory','directory for mass flux profiles');
			fielddisplay(self,'mass_flux_profiles','list of mass_flux profiles');
			fielddisplay(self,'mass_flux_segments','');
			fielddisplay(self,'adjacency','');
			fielddisplay(self,'vertex_weight','weight applied to each mesh vertex');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','isdakota','format','Boolean');
			if ~self.isdakota, 
				WriteData(fid,prefix,'data',false,'name','md.qmu.mass_flux_segments_present','format','Boolean');
				return; 
			end
			WriteData(fid,prefix,'object',self,'fieldname','partition','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'fieldname','numberofpartitions','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','numberofresponses','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','variabledescriptors','format','StringArray');
			WriteData(fid,prefix,'object',self,'fieldname','responsedescriptors','format','StringArray');
			if ~isempty(self.mass_flux_segments), 
				WriteData(fid,prefix,'data',self.mass_flux_segments,'name','md.qmu.mass_flux_segments','format','MatArray');
				flag=true; 
			else 
				flag=false; 
			end
			WriteData(fid,prefix,'data',flag,'name','md.qmu.mass_flux_segments_present','format','Boolean');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			if self.isdakota,
				error('qmu savemodeljs error message: not supported yet!');
			end

		end % }}}
	end
end
