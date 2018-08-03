function ispresent=waitonlock(md)
%WAITONLOCK - wait for a file
%
%   This routine will return when a file named 'lockfilename' is written to disk.
%   Also check for outlog file be cause it might bewritten several seconds
%   after the lock file.
%   If the time limit given in input is exceeded, return 0
%
%   Usage:
%      flag=waitonlock(md)

%Return if waitonlock < 0 (no need to wait)

%Get lockfilename (lock file) and options
executionpath = md.cluster.executionpath;
timelimit     = md.settings.waitonlock;
cluster       = md.cluster;

if isa(cluster,'pfe') && cluster.interactive>0
	lockfilename  = [executionpath '/Interactive' num2str(cluster.interactive) '/' md.miscellaneous.name '.lock'];
	logfilename   = [executionpath '/Interactive' num2str(cluster.interactive) '/' md.miscellaneous.name '.outlog'];
else
	lockfilename  = [executionpath '/' md.private.runtimename '/' md.miscellaneous.name '.lock'];
	logfilename   = [executionpath '/' md.private.runtimename '/' md.miscellaneous.name '.outlog'];
end


%If we are using the generic cluster in interactive mode, job is already complete
if (isa(cluster,'generic') & cluster.interactive) | isa(cluster,'generic_static'),
	%We are in interactive mode, no need to check for job completion
	ispresent=1;
	return;
end

%initialize time and file presence test flag
time=0; ispresent=0; time0=clock;
disp(['waiting for ' lockfilename ' hold on... (Ctrl+C to exit)'])

%prepare command if the job is not running on the local machine
if ~strcmpi(oshostname(),cluster.name),
	login = cluster.login;
	port  = 0;
	if isprop(cluster,'port') 
		port = cluster.port;
	end
	if port,
		command = ['ssh -l ' login ' -p ' num2str(port) ' localhost "[ -f ' lockfilename ' ] && [ -f ' logfilename ' ]" 2>/dev/null'];
	elseif isa(cluster,'cloud')
		command = [' [ -f ' lockfilename ' ] && [ -f ' logfilename ' ] 2>/dev/null'];
		command = [starcluster() ' sshmaster ' cluster.name ' --user ' cluster.login ' ''' command ''''];
	else
		command = ['ssh -l ' login ' ' cluster.name ' "[ -f ' lockfilename ' ] && [ -f ' logfilename ' ]" 2>/dev/null'];
	end
end

%loop till file .lock exist or time is up
while (ispresent==0 & time<timelimit)
	if strcmpi(oshostname(),cluster.name),
		pause(1);
		ispresent=(exist(lockfilename,'file') & exist(logfilename,'file'));
		time=etime(clock,time0)/60;
	else
		pause(5);
		time=etime(clock,time0);
		fprintf('\rchecking for job completion (time: %i min %i sec)      ',floor(time/60),floor(rem(time,60)));
		time=time/60; %converts time from sec to min
		ispresent=~system(command);
		if ispresent, fprintf('\n'); end
	end
end

%build output
if (time>timelimit),
	disp('Time limit exceeded. Increase md.settings.waitonlock');
	disp('The results must be loaded manually with md=loadresultsfromcluster(md).');
	error(['waitonlock error message: time limit exceeded']);
end
