%GENERIC cluster class definition
%
%   Usage:
%      cluster=generic_static('name','astrid','np',3);

classdef generic_static
	properties (SetAccess=public) 
		% {{{
		name='';
		np=1;
		codepath=fileparts(which('issm.exe'));
		executionpath = '.';
		interactive = 1;
		shell='/bin/sh';
		%}}}
	end
	methods
		function cluster=generic_static(varargin) % {{{

			%use provided options to change fields
			options=pairoptions(varargin{:});

			%get name
			cluster.name=getfieldvalue(options,'name',oshostname());

			%initialize cluster using user settings if provided
			if (exist([cluster.name '_settings'])==2), eval([cluster.name '_settings']); end

			%OK get other fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
		function disp(cluster) % {{{
			%  display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name: %s',cluster.name));
			disp(sprintf('    np: %i',cluster.np));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    shell: %s',cluster.shell));
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{
			if cluster.np<1
				md = checkmessage(md,['number of processors should be at least 1']);
			end
			if isnan(cluster.np),
				md = checkmessage(md,'number of processors should not be NaN!');
			end
		end
		%}}}
		function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			%Check that issm.exe exists in the right path
			if ~exist([cluster.codepath '/issm.exe'],'file'),
				error(['File ' cluster.codepath '/issm.exe does not exist']);
			end

			%Now process codepath and replace empty spaces with \ to avoid errors in queuing script
			codepath2=strrep(cluster.codepath,' ','\ ');

			%write queuing script
			%what is the executable being called?
			executable='issm.exe';
			if isdakota,
				version=IssmConfig('_DAKOTA_VERSION_'); version=str2num(version(1:3));
				if (version>=6),
					executable='issm_dakota.exe';
				end
			end
			if isoceancoupling,
				executable='issm_ocean.exe';
			end

			%write queuing script 
			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#!%s\n',cluster.shell);
			fprintf(fid,['%s/mpiexec -np %i %s/%s %s %s %s \n'],codepath2,cluster.np,codepath2,executable,solution,'./',modelname);
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			fid=fopen([modelname '.errlog'],'w'); fclose(fid);
			fid=fopen([modelname '.outlog'],'w'); fclose(fid);
		end
		%}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist)% {{{

			%do nothing
		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch)% {{{

			if ~ispc,

				%figure out what shell extension we will use:
				if isempty(strfind(cluster.shell,'csh')),
					shellext='sh';
				else
					shellext='csh';
				end

				disp('launching solution sequence');
				launchcommand=['source  ' modelname '.queue '];
				issmssh(cluster.name,'',0,launchcommand);
			else
				system([modelname '.bat']);
			end
		end %}}}
		function Download(cluster,dirname,filelist)% {{{
				%do nothing
				return;
		end %}}}
	end
end
