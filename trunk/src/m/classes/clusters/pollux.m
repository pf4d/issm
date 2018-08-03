%POLLUX class definition
%
%   Usage:
%      cluster=pollux();
%      cluster=pollux('np',3);
%      cluster=pollux('np',3,'login','username');

classdef pollux
    properties (SetAccess=public) 
		 % {{{
		 name='pollux'
		 login='username';
		 np=128;
		 port=0;
		 queue='shortp';
		 time=180;
		 codepath='/workc/edw/issm-2.0/bin'
		 executionpath='/workc/edw/Testing/Execution'
		 %}}}
	 end
	 methods
		 function cluster=pollux(varargin) % {{{
			 cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		 end
		 %}}}
		 function disp(cluster) % {{{
			 %  display the object
			 disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			 disp(sprintf('    name: %s',cluster.name));
			 disp(sprintf('    login: %s',cluster.login));
			 disp(sprintf('    port: %i',cluster.port));
			 disp(sprintf('    np: %i',cluster.np));
			 disp(sprintf('    queue: %s',cluster.queue));
			 disp(sprintf('    time: %i',cluster.time));
			 disp(sprintf('    codepath: %s',cluster.codepath));
			 disp(sprintf('    executionpath: %s',cluster.executionpath));
		 end
		 %}}}
		 function md = checkconsistency(cluster,md,solution,analyses) % {{{

			 available_queues={'shortp','longp'};
			 queue_requirements_time=[180 720];
			 queue_requirements_np=[128 128];

			 QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.np,cluster.time)
		 end
		 %}}}
		 function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			 if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			 if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#!/bin/sh\n');
			 fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); %walltime is in seconds.
			 fprintf(fid,'#PBS -N %s\n',modelname);
			 fprintf(fid,'#PBS -l ncpus=%i\n',cluster.np);
			 if ~isempty(queue),
				 fprintf(fid,'#PBS -q %s\n',cluster.queue);
			 end
			 fprintf(fid,'#PBS -o %s.outlog \n',modelname);
			 fprintf(fid,'#PBS -e %s.errlog \n',modelname);
			 fprintf(fid,'export PBS_O_WORKDIR=%s\n',[cluster.executionpath '/' dirname]);
			 fprintf(fid,'cd $PBS_O_WORKDIR\n');
			 fprintf(fid,'export OMP_NUM_THREADS=1\n');
			 fprintf(fid,'dplace -s1 -c0-%i mpiexec -np %i %s/issm.exe %s %s %s',cluster.np-1,cluster.np,cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
			 fclose(fid);

		 end
		 %}}}
		 function LaunchQueueJob(cluster,modelname,dirname,filelist,restart)% {{{

			 disp('launching solution sequence on remote cluster');
			 if ~isempty(restart)
				 launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && qsub ' modelname '.queue '];
			 else
				 launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					 ' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && qsub ' modelname '.queue '];
			 end
			 issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		 end %}}}
		 function UploadQueueJob(cluster,modelname,dirname,filelist)% {{{

			 %compress the files into one zip.
			 compressstring=['tar -zcf ' dirname '.tar.gz '];
			 for i=1:numel(filelist),
				 compressstring = [compressstring ' ' filelist{i}];
			 end
			 if cluster.interactive,
				 compressstring = [compressstring ' ' modelname '.errlog ' modelname '.outlog '];
			 end
			 system(compressstring);

		 end %}}}
		 function Download(cluster,dirname,filelist)% {{{

			 %copy files from cluster to current directory
			 directory=[cluster.executionpath '/' dirname '/'];
			 issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);

		 end %}}}
	end
end
