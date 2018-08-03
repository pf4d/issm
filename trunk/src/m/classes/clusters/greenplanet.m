%PFE class definition
%
%   Usage:
%      cluster=greenplanet();
%      cluster=greenplanet('np',3);
%      cluster=greenplanet('np',3,'login','username');

classdef greenplanet
    properties (SetAccess=public)  
		 % {{{
		 name          = 'greenplanet'
		 login         = '';
		 numnodes      = 20;
		 cpuspernode   = 8;
		 port          = 8000;
		 queue         = 'c6145';
		 codepath      = '';
		 executionpath = '';
		 interactive   = 0;
		 time          = 24*60;
		 memory        = 2;
	 end
	 %}}}
	 methods
		 function cluster=greenplanet(varargin) % {{{

			 %initialize cluster using default settings if provided
			 if (exist('greenplanet_settings')==2), greenplanet_settings; end

			 %use provided options to change fields
			 cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		 end
		 %}}}
		 function disp(cluster) % {{{
			 %  display the object
			 disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			 disp(sprintf('    name: %s',cluster.name));
			 disp(sprintf('    login: %s',cluster.login));
			 disp(sprintf('    port: %i',cluster.port));
			 disp(sprintf('    numnodes: %i',cluster.numnodes));
			 disp(sprintf('    cpuspernode: %i',cluster.cpuspernode));
			 disp(sprintf('    np: %i',cluster.cpuspernode*cluster.numnodes));
			 disp(sprintf('    queue: %s',cluster.queue));
			 disp(sprintf('    codepath: %s',cluster.codepath));
			 disp(sprintf('    executionpath: %s',cluster.executionpath));
			 disp(sprintf('    interactive: %i',cluster.interactive));
			 disp(sprintf('    time: %i',cluster.time));
			 disp(sprintf('    memory: %i',cluster.memory));
		 end
		 %}}}
		 function numprocs=np(cluster) % {{{
			 %compute number of processors
			 numprocs=cluster.numnodes*cluster.cpuspernode;
		 end
		 %}}}
		 function md = checkconsistency(cluster,md,solution,analyses) % {{{

			 available_queues={'c6145','default'};
			 queue_requirements_time=[Inf Inf];
			 queue_requirements_np=[80 80];

			 QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.np,1)

			 %Miscelaneous
			 if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			 if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			 if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end

		 end
		 %}}}
		 function BuildKrigingQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			 if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			 if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#!/bin/bash\n');
			 fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
			 fprintf(fid,'#SBATCH -p %s \n',cluster.queue);
			 fprintf(fid,'#SBATCH -N %i -n %i\n',cluster.numnodes,cluster.cpuspernode);
			 fprintf(fid,'#SBATCH --time=%i\n',cluster.time*60); %walltime is in seconds.
			 fprintf(fid,'#SBATCH --mem-per-cpu=%igb\n',cluster.memory);
			 fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			 fprintf(fid,'#SBATCH -e %s.errlog \n\n',modelname);
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			 fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			 fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			 fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s\n',cluster.np,cluster.codepath,[cluster.executionpath '/' modelname],modelname);
			 if ~io_gather, %concatenate the output files:
				 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			 end
			 fclose(fid);
		 end
		 %}}}
		 function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota) % {{{

			 if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			 if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#!/bin/bash\n');
			 fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
			 fprintf(fid,'#SBATCH -p %s \n',cluster.queue);
			 fprintf(fid,'#SBATCH -N %i -n %i\n',cluster.numnodes,cluster.cpuspernode);
			 fprintf(fid,'#SBATCH --time=%i\n',cluster.time*60); %walltime is in seconds.
			 fprintf(fid,'#SBATCH --mem-per-cpu=%igb\n',cluster.memory);
			 fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			 fprintf(fid,'#SBATCH -e %s.errlog \n\n',modelname);
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			 fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			 fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			 fprintf(fid,'mpiexec -np %i %s/issm.exe %s %s %s\n',cluster.np,cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
			 if ~io_gather, %concatenate the output files:
				 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			 end
			 fclose(fid);

			 %in interactive mode, create a run file, and errlog and outlog file
			 if cluster.interactive,
				 fid=fopen([modelname '.run'],'w');
				 fprintf(fid,'mpiexec -np %i %s/issm.exe %s %s %s\n',cluster.np,cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
				 if ~io_gather, %concatenate the output files:
					 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
				 end
				 fclose(fid);
				 fid=fopen([modelname '.errlog'],'w');
				 fclose(fid);
				 fid=fopen([modelname '.outlog'],'w');
				 fclose(fid);
			 end
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

			 disp('uploading input file and queueing script');
			 issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']});

		 end %}}}
		 function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch)% {{{

			 disp('launching solution sequence on remote cluster');
			 if ~isempty(restart)
				 launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && qsub ' modelname '.queue '];
			 else
				 launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					 ' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && qsub ' modelname '.queue '];
			 end
			 issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		 end %}}}
		 function Download(cluster,dirname,filelist)% {{{

			 %copy files from cluster to current directory
			 directory=[cluster.executionpath '/' dirname '/'];
			 issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);

		 end %}}}
	end
end
