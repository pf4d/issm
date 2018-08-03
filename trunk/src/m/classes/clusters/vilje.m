%PFE class definition
%
%   Usage:
%      cluster=greenplanet();
%      cluster=greenplanet('np',3);
%      cluster=greenplanet('np',3,'login','username');

classdef vilje
    properties (SetAccess=public)  
		 % {{{
     name           = 'vilje';
		 login          = '';
		 numnodes       = 2;
		 cpuspernode    = 32;
     procspernodes  = 16;
     mem            = 28;
		 queue          = 'workq';
		 time           = 2*60;
		 codepath       = '';
		 executionpath  = '';
		 interactive    = 0;
     port           = [];
     accountname    = '';

	 end
	 properties (SetAccess=private) 
		 np=numnodes*procspernodes;
		 % }}}
	 end
	 methods
		 function cluster=vilje(varargin) % {{{

			 %initialize cluster using default settings if provided
			 if (exist('vilje_settings')==2), vilje_settings; end

			 %use provided options to change fields
			 cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		 end
		 %}}}
		 function disp(cluster) % {{{
			 %  display the object
			 disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			 disp(sprintf('    name: %s',cluster.name));
			 disp(sprintf('    login: %s',cluster.login));
       disp(sprintf('    accountname: %s',cluster.accountname));
			 disp(sprintf('    numnodes: %i',cluster.numnodes));
			 disp(sprintf('    cpuspernode: %i',cluster.cpuspernode));
			 disp(sprintf('    np: %i', cluster.cpuspernode*cluster.numnodes));
			 disp(sprintf('    procspernodes: %i',cluster.procspernodes));
			 disp(sprintf('    queue: %s',cluster.queue));
			 disp(sprintf('    codepath: %s',cluster.codepath));
			 disp(sprintf('    executionpath: %s',cluster.executionpath));
			 disp(sprintf('    interactive: %i',cluster.interactive));
			 disp(sprintf('    time: %i',cluster.time));
			 disp(sprintf('    memory: %i',cluster.mem));
		 end
		 %}}}
		 function md = checkconsistency(cluster,md,solution,analyses) % {{{

			 available_queues={'workq'};
			 queue_requirements_time=[5*24*60];
			 queue_requirements_np=[30];

			 QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.np,1)

			 %Miscelaneous
			 if isempty(cluster.login), md = checkmessage(md,'login empty'); end
       if isempty(cluster.accountname), md = checkmessage(md,'accountname empty'); end
			 if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			 if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end

		 end
		 %}}}
		 function BuildKrigingQueueScript(cluster,modelname,solution,io_gather,isvalgrind,isgprof) % {{{

			 if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			 if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			 %compute number of processors
			 cluster.np=cluster.numnodes*cluster.cpuspernode;

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#PBS -S /bin/bash\n');
			 fprintf(fid,'#PBS -N %s\n',modelname);
       fprintf(fid,'#PBS -l select=%i:ncpus=%i:mpiprocs=%i\n',cluster.numnodes,cluster.cpuspernode,16);
       fprintf(fid,'#PBS -l walltime=%s\n',cluster.time); %walltime is in seconds.
       fprintf(fid,'#PBS -A %s\n',cluster.accountname);
			 fprintf(fid,'#PBS -o %s.outlog \n',modelname);
			 fprintf(fid,'#PBS -e %s.errlog \n\n',modelname);
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			 fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			 fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,modelname);
			 fprintf(fid,'mpiexec_mpt -n %i %s/kriging.exe %s %s\n',cluster.np,cluster.codepath,[cluster.executionpath '/' modelname],modelname);
			 if ~io_gather, %concatenate the output files:
				 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			 end
			 fclose(fid);
		 end
		 %}}}
		 function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			 if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			 if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

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

			 %compute number of processors
			 cluster.np=cluster.numnodes*cluster.cpuspernode;                     
       shortname = substring(modelname,1,min(12,length(modelname)));

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#PBS -S /bin/bash\n');
			 fprintf(fid,'#PBS -N %s\n',shortname);
			 fprintf(fid,'#PBS -q %s \n',cluster.queue);
       fprintf(fid,'#PBS -l select=%i:ncpus=%i:mpiprocs=%i\n',cluster.numnodes,cluster.cpuspernode,cluster.procspernodes);
			 fprintf(fid,'#PBS -l walltime=%s\n',duration(0,cluster.time,0)); %walltime is in minutes.
       fprintf(fid,'#PBS -A %s\n',cluster.accountname);
			 fprintf(fid,'#PBS -o %s.outlog \n',[cluster.executionpath '/' dirname '/' modelname]);
			 fprintf(fid,'#PBS -e %s.errlog \n\n',[cluster.executionpath '/' dirname '/' modelname]);
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			 fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			 fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
       fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);

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
			 system(compressstring);
			 disp('uploading input file and queueing script');
			 directory=cluster.executionpath;
			 issmbbftpout(cluster.name,directory,cluster.login,cluster.port,cluster.numstreams,{[dirname '.tar.gz']});

		 end
		 %}}}


		 function LaunchQueueJob(cluster,modelname,dirname,filelist)% {{{

			 disp('launching solution sequence on remote cluster');
			  launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
											 ' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && qsub ' modelname '.queue '];
			 issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		 end %}}}

		 function Download(cluster,dirname,filelist)% {{{

			 %copy files from cluster to current directory
			 directory=[cluster.executionpath '/' dirname '/'];
			 issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);

		 end %}}}
	end
end
