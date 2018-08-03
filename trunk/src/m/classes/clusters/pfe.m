%PFE class definition
%
%   Usage:
%      cluster=pfe();
%      cluster=pfe('np',3);
%      cluster=pfe('np',3,'login','username');

classdef pfe
    properties (SetAccess=public)  
		 % {{{
		 name           = 'pfe'
		 login          = '';
		 modules        = {'comp-intel/2016.2.181' 'mpi-sgi/mpt'};
		 numnodes       = 20;
		 cpuspernode    = 8;
		 port           = 1025;
		 queue          = 'long';
		 time           = 12*60;
		 processor      = 'ivy';
		 codepath       = '';
		 executionpath  = '';
		 grouplist     = 's1690';
		 interactive    = 0;
		 bbftp          = 0;
		 numstreams     = 8;
		 hyperthreading = 0;
	 end
	 %}}}
	 methods
		 function cluster=pfe(varargin) % {{{

			 %initialize cluster using default settings if provided
			 if (exist('pfe_settings')==2), pfe_settings; end

			 %use provided options to change fields
			 cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		 end
		 %}}}
		 function disp(cluster) % {{{
			 %  display the object
			 disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			 disp(sprintf('    name: %s',cluster.name));
			 disp(sprintf('    login: %s',cluster.login));
			 modules=''; for i=1:length(cluster.modules), modules=[modules cluster.modules{i} ',']; end; modules=modules(1:end-1); 
			 disp(sprintf('    modules: %s',modules));
			 disp(sprintf('    port: %i',cluster.port));
			 disp(sprintf('    numnodes: %i',cluster.numnodes));
			 disp(sprintf('    cpuspernode: %i',cluster.cpuspernode));
			 disp(sprintf('    np: %i',cluster.cpuspernode*cluster.numnodes));
			 disp(sprintf('    queue: %s',cluster.queue));
			 disp(sprintf('    time: %i',cluster.time));
			 disp(sprintf('    processor: %s',cluster.processor));
			 disp(sprintf('    codepath: %s',cluster.codepath));
			 disp(sprintf('    executionpath: %s',cluster.executionpath));
			 disp(sprintf('    grouplist: %s',cluster.grouplist));
			 disp(sprintf('    interactive: %i',cluster.interactive));
			 disp(sprintf('    hyperthreading: %i',cluster.hyperthreading));
		 end
		 %}}}
		 function numprocs=np(cluster) % {{{
			 %compute number of processors
			 numprocs=cluster.numnodes*cluster.cpuspernode;
		 end
		 %}}}
		 function md = checkconsistency(cluster,md,solution,analyses) % {{{

			 available_queues={'long','normal','debug','devel','alphatst@pbspl233'};
			 queue_requirements_time=[5*24*60 8*60 2*60 2*60 24*60];
			 queue_requirements_np=[2048 2048 150 150 2048];

			 QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.np,cluster.time)

			 %now, check cluster.cpuspernode according to processor type
			 if strcmpi(cluster.processor,'wes'),
				 if cluster.hyperthreading,
					 if ((cluster.cpuspernode>24 ) | (cluster.cpuspernode<1)),
						 md = checkmessage(md,'cpuspernode should be between 1 and 24 for ''wes'' processors in hyperthreading mode');
					 end
				 else
					 if ((cluster.cpuspernode>12 ) | (cluster.cpuspernode<1)),
						 md = checkmessage(md,'cpuspernode should be between 1 and 12 for ''wes'' processors');
					 end
				 end
			 elseif strcmpi(cluster.processor,'ivy'),
				 if cluster.hyperthreading,
					 if ((cluster.cpuspernode>40 ) | (cluster.cpuspernode<1)),
						 md = checkmessage(md,'cpuspernode should be between 1 and 40 for ''ivy'' processors in hyperthreading mode');
					 end
				 else
					 if ((cluster.cpuspernode>20 ) | (cluster.cpuspernode<1)),
						 md = checkmessage(md,'cpuspernode should be between 1 and 20 for ''ivy'' processors');
					 end
				 end
			 elseif strcmpi(cluster.processor,'bro'),
				 if cluster.hyperthreading,
					 if ((cluster.cpuspernode>56 ) | (cluster.cpuspernode<1)),
						 md = checkmessage(md,'cpuspernode should be between 1 and 56 for ''bro'' processors in hyperthreading mode');
					 end
				 else
					 if ((cluster.cpuspernode>28 ) | (cluster.cpuspernode<1)),
						 md = checkmessage(md,'cpuspernode should be between 1 and 28 for ''bro'' processors');
					 end
				 end
			 elseif strcmpi(cluster.processor,'has'),
				 if cluster.hyperthreading,
					 if ((cluster.cpuspernode>48 ) | (cluster.cpuspernode<1)),
						 md = checkmessage(md,'cpuspernode should be between 1 and 48 for ''has'' processors in hyperthreading mode');
					 end
				 else
					 if ((cluster.cpuspernode>24 ) | (cluster.cpuspernode<1)),
						 md = checkmessage(md,'cpuspernode should be between 1 and 24 for ''has'' processors');
					 end
				 end
			 
			 elseif strcmpi(cluster.processor,'san'),
				 if cluster.hyperthreading,
					 if ((cluster.cpuspernode>32 ) | (cluster.cpuspernode<1)),
						 md = checkmessage(md,'cpuspernode should be between 1 and 32 for ''san'' processors in hyperthreading mode');
					 end
				 else
					 if ((cluster.cpuspernode>16 ) | (cluster.cpuspernode<1)),
						 md = checkmessage(md,'cpuspernode should be between 1 and 16 for ''san'' processors');
					 end
				 end

			 else
				 md = checkmessage(md,'unknown processor type, should be ''wes'' or ''has'' or ''ivy'' or ''san''');
			 end

			 %Miscelaneous
			 if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			 if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			 if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end
			 if isempty(cluster.grouplist), md = checkmessage(md,'grouplist empty'); end

		 end
		 %}}}
		 function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

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

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#PBS -S /bin/bash\n');
%			 fprintf(fid,'#PBS -N %s\n',modelname);
			 fprintf(fid,'#PBS -l select=%i:ncpus=%i:model=%s\n',cluster.numnodes,cluster.cpuspernode,cluster.processor);
			 fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); %walltime is in seconds.
			 fprintf(fid,'#PBS -q %s \n',cluster.queue);
			 fprintf(fid,'#PBS -W group_list=%s\n',cluster.grouplist);
			 fprintf(fid,'#PBS -m e\n');
			 fprintf(fid,'#PBS -o %s.outlog \n',[cluster.executionpath '/' dirname '/' modelname]);
			 fprintf(fid,'#PBS -e %s.errlog \n\n',[cluster.executionpath '/' dirname '/' modelname]);
			 fprintf(fid,'. /usr/share/modules/init/bash\n\n');
			 fprintf(fid,'module load comp-intel/2016.2.181\n');
			 fprintf(fid,'module load mpi-sgi/mpt\n');
			 fprintf(fid,'export PATH="$PATH:."\n\n');
			 fprintf(fid,'export MPI_GROUP_MAX=64\n\n');
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			 fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			 fprintf(fid,'cd %s/%s/\n\n',cluster.executionpath,dirname);
			 if ~isvalgrind,
				 fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			 else
				 fprintf(fid,'mpiexec -np %i valgrind --leak-check=full %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			 end
			 if ~io_gather, %concatenate the output files:
				 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			 end
			 fclose(fid);

			 %in interactive mode, create a run file, and errlog and outlog file
			 if cluster.interactive,
				 fid=fopen([modelname '.run'],'w');
				 if cluster.interactive==10,
						 fprintf(fid,'module unload mpi-mvapich2/1.4.1/gcc\n');
						 fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[pwd() '/run'],modelname);
				 else
					 if ~isvalgrind,
						 fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[cluster.executionpath '/Interactive' num2str(cluster.interactive)],modelname);
					 else
						 fprintf(fid,'mpiexec -np %i valgrind --leak-check=full %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[cluster.executionpath '/Interactive' num2str(cluster.interactive)],modelname);
					 end
				 end
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
		 function BuildQueueScriptMultipleModels(cluster,dirname,modelname,solution,dirnames,modelnames,nps) % {{{

			 %some checks: 
			 if isempty(modelname), error('BuildQueueScriptMultipleModels error message: need a non empty model name!');end

			 %what is the executable being called? 
			 executable='issm_slr.exe';

			 if ispc(), error('BuildQueueScriptMultipleModels not support yet on windows machines');end;

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');

			 fprintf(fid,'#PBS -S /bin/bash\n');
			 fprintf(fid,'#PBS -N %s\n',modelname);
			 fprintf(fid,'#PBS -l select=%i:ncpus=%i:model=%s\n',cluster.numnodes,cluster.cpuspernode,cluster.processor);
			 fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); %walltime is in seconds.
			 fprintf(fid,'#PBS -q %s \n',cluster.queue);
			 fprintf(fid,'#PBS -W group_list=%s\n',cluster.grouplist);
			 fprintf(fid,'#PBS -m e\n');
			 fprintf(fid,'#PBS -o %s.outlog \n',[cluster.executionpath '/' dirname '/' modelname]);
			 fprintf(fid,'#PBS -e %s.errlog \n\n',[cluster.executionpath '/' dirname '/' modelname]);
			 fprintf(fid,'. /usr/share/modules/init/bash\n\n');
			 fprintf(fid,'module load comp-intel/2016.2.181\n');
			 fprintf(fid,'module load mpi-sgi/mpt\n');
			 fprintf(fid,'export PATH="$PATH:."\n\n');
			 fprintf(fid,'export MPI_GROUP_MAX=64\n\n');
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			 fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			 fprintf(fid,'cd %s/%s/\n\n',cluster.executionpath,dirname);

			 %number of cpus: 
			 mpistring=sprintf('mpiexec -np %i ',cluster.numnodes*cluster.cpuspernode);

			 %executable: 
			 mpistring=[mpistring sprintf('%s/%s ',cluster.codepath,executable)];

			 %solution name: 
			 mpistring=[mpistring sprintf('%s ',solution)];

			 %execution directory and model name: 
			 mpistring=[mpistring sprintf('%s/%s %s',cluster.executionpath,dirname,modelname)];

			 %inform main executable of how many icecaps, glaciers and earth models are being run: 
			 mpistring=[mpistring sprintf(' %i ',length(dirnames))];

			 %icecaps, glaciers and earth location, names and number of processors associated:
			 for i=1:length(dirnames),
				 mpistring=[mpistring sprintf(' %s/%s %s %i ',cluster.executionpath,dirnames{i},modelnames{i},nps{i})];
			 end

			 %write this long string to disk: 
			 fprintf(fid,mpistring);
			 fclose(fid);

			 if cluster.interactive,
				 fid=fopen([modelname '.run'],'w');

				 %number of cpus: 
				 mpistring=sprintf('mpiexec -np %i ',cluster.numnodes*cluster.cpuspernode);

				 %executable: 
				 mpistring=[mpistring sprintf('%s/%s ',cluster.codepath,executable)];

				 %solution name: 
				 mpistring=[mpistring sprintf('%s ',solution)];

				 %execution directory and model name: 
				 mpistring=[mpistring sprintf('%s/%s %s',cluster.executionpath,dirname,modelname)];

				 %inform main executable of how many icecaps, glaciers and earth models are being run: 
				 mpistring=[mpistring sprintf(' %i ',length(dirnames))];

				 %icecaps, glaciers and earth location, names and number of processors associated:
				 for i=1:length(dirnames),
					 mpistring=[mpistring sprintf(' %s/Interactive%i %s %i ',cluster.executionpath,cluster.interactive,modelnames{i},nps{i})];
				 end

				 %write this long string to disk: 
				 fprintf(fid,mpistring);
				 fclose(fid);

				 fid=fopen([modelname '.errlog'],'w');
				 fclose(fid);
				 fid=fopen([modelname '.outlog'],'w');
				 fclose(fid);
			 end
		 end
		 %}}}
		 function BuildKrigingQueueScript(cluster,modelname,solution,io_gather,isvalgrind,isgprof) % {{{

			 if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#PBS -S /bin/bash\n');
			 %			 fprintf(fid,'#PBS -N %s\n',modelname);
			 fprintf(fid,'#PBS -l select=%i:ncpus=%i:model=%s\n',cluster.numnodes,cluster.cpuspernode,cluster.processor);
			 fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); %walltime is in seconds.
			 fprintf(fid,'#PBS -q %s \n',cluster.queue);
			 fprintf(fid,'#PBS -W group_list=%s\n',cluster.grouplist);
			 fprintf(fid,'#PBS -m e\n');
			 fprintf(fid,'#PBS -o %s.outlog \n',modelname);
			 fprintf(fid,'#PBS -e %s.errlog \n\n',modelname);
			 fprintf(fid,'. /usr/share/modules/init/bash\n\n');
			 for i=1:numel(cluster.modules),
				 fprintf(fid,['module load ' cluster.modules{i} '\n']);
			 end
			 fprintf(fid,'export PATH="$PATH:."\n');
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			 fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			 fprintf(fid,'export MPI_GROUP_MAX=64\n\n');
			 fprintf(fid,'cd %s/%s/\n\n',cluster.executionpath,modelname);
			 fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s\n',cluster.np,cluster.codepath,[cluster.executionpath '/' modelname],modelname); %FIXME
			 if ~io_gather, %concatenate the output files:
				 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			 end
			 fclose(fid);

			 %in interactive mode, create a run file, and errlog and outlog file
			 if cluster.interactive,
				 fid=fopen([modelname '.run'],'w');
				 if ~isvalgrind,
					 fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s\n',cluster.np,cluster.codepath,[cluster.executionpath '/' modelname],modelname);
				 else
					 fprintf(fid,'mpiexec -np %i valgrind --leak-check=full %s/kriging.exe %s %s\n',cluster.np,cluster.codepath,[cluster.executionpath '/' modelname],modelname);
				 end
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
		 function BuildOceanQueueScript(np,cluster,modelname) % {{{

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#PBS -S /bin/bash\n');
			 fprintf(fid,'#PBS -l select=1:ncpus=%i:model=%s\n',np,cluster.processor);
			 fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); %walltime is in seconds.
			 fprintf(fid,'#PBS -q %s \n',cluster.queue);
			 fprintf(fid,'#PBS -W group_list=%s\n',cluster.grouplist);
			 fprintf(fid,'#PBS -m e\n');
			 fprintf(fid,'#PBS -o %s.outlog \n',modelname);
			 fprintf(fid,'#PBS -e %s.errlog \n\n',modelname);
			 fprintf(fid,'. /usr/share/modules/init/bash\n\n');
			 fprintf(fid,'module load comp-intel/2015.0.090\n');
			 fprintf(fid,'module load test/mpt.2.11r8\n');
			 fprintf(fid,'module load netcdf/4.0\n');
			 fprintf(fid,'module load mpi-mvapich2/1.4.1/gcc\n');
			 fprintf(fid,'module load gcc/4.4.4\n');
			 fprintf(fid,'export PATH="$PATH:."\n');
			 fprintf(fid,'export MPI_GROUP_MAX=64\n\n');
			 fprintf(fid,['cd ' pwd() ' \n\n']);
			 fprintf(fid,'mpiexec -np %i ./mitgcmuv\n',np); 
		%	 if ~io_gather, %concatenate the output files:
		%		 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
		%	 end
			 fclose(fid);

			 %in interactive mode, create a run file, and errlog and outlog file
			 if cluster.interactive,
				 fid=fopen([modelname '.run'],'w');
				 fprintf(fid,'module load mpi-mvapich2/1.4.1/gcc\n');
				 fprintf(fid,['mpiexec -np %i ./mitgcmuv \n'],np);
				 fprintf(fid,['touch ' modelname '.lock %s\n']);
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
				 compressstring = [compressstring ' ' modelname '.run '  modelname '.errlog ' modelname '.outlog '];
			 end
			 system(compressstring);

			 disp('uploading input file and queueing script');
			 if cluster.interactive==10,
				 directory=[pwd() '/run/'];
			 elseif cluster.interactive,
				 directory=[cluster.executionpath '/Interactive' num2str(cluster.interactive)];
			 else 
				 directory=cluster.executionpath;
			 end

			 if ~cluster.bbftp,
				 issmscpout(cluster.name,directory,cluster.login,cluster.port,{[dirname '.tar.gz']});
			 else
				 issmbbftpout(cluster.name,directory,cluster.login,cluster.port,cluster.numstreams,{[dirname '.tar.gz']});
			 end

		 end
		 %}}}
		 function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch)% {{{

			 %lauch command, to be executed via ssh
			 if ~cluster.interactive, 
				 if ~isempty(restart)
					 launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && qsub ' modelname '.queue '];
				 else
					 launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						 ' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && qsub ' modelname '.queue '];
				 end
			 else
				 if ~isempty(restart)
					 launchcommand=['cd ' cluster.executionpath '/Interactive' num2str(cluster.interactive)];
				 else
					 if cluster.interactive==10,
						 launchcommand=['cd ' pwd() '/run && tar -zxf ' dirname '.tar.gz'];
					 else
						 launchcommand=['cd ' cluster.executionpath '/Interactive' num2str(cluster.interactive) ' && tar -zxf ' dirname '.tar.gz'];
					 end
				 end
			 end

			 disp('launching solution sequence on remote cluster');
			 issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		 end
		 %}}}
		 function Download(cluster,dirname,filelist)% {{{

			 %copy files from cluster to current directory
			 if cluster.interactive==10,
				 directory=[pwd() '/run/'];
			 elseif ~cluster.interactive,
				 directory=[cluster.executionpath '/' dirname '/'];
			 else
				 directory=[cluster.executionpath '/Interactive' num2str(cluster.interactive) '/'];
			 end

			 if ~cluster.bbftp,
				 issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);
			 else
				 issmbbftpin(cluster.name, cluster.login, cluster.port, cluster.numstreams, directory, filelist);
			 end

		 end %}}}
	end
end
