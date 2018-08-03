function md=loadresultsfromcluster(md,runtimename)
%LOADRESULTSFROMCLUSTER - load results of solution sequence from cluster
%
%   Usage:
%      md=loadresultsfromcluster(md,runtimename);

%retrieve cluster, to be able to call its methods
cluster=md.cluster;

if nargin==2,
	md.private.runtimename=runtimename;
end

%Download outputs from the cluster
filelist={[md.miscellaneous.name '.outlog'],[md.miscellaneous.name '.errlog']};
if md.qmu.isdakota,
	filelist{end+1}=[md.miscellaneous.name '.qmu.err'];
	filelist{end+1}=[md.miscellaneous.name '.qmu.out'];
	if isfield(md.qmu.params,'tabular_graphics_data'),
		if md.qmu.params.tabular_graphics_data==true,
			filelist{end+1}='dakota_tabular.dat';
		end
	end
else
	filelist{end+1}=[md.miscellaneous.name '.outbin'];
end
Download(cluster,md.private.runtimename,filelist);

%If we are here, no errors in the solution sequence, call loadresultsfromdisk.
md=loadresultsfromdisk(md,[md.miscellaneous.name '.outbin']);

%erase the log and output files
if md.qmu.isdakota,
	delete([['qmu' num2str(feature('GetPid')) '/'] md.miscellaneous.name '.outlog']);
	delete([['qmu' num2str(feature('GetPid')) '/']  md.miscellaneous.name '.errlog']);
else
	delete([md.miscellaneous.name '.outlog']);
	delete([md.miscellaneous.name '.errlog']);
	delete([md.miscellaneous.name '.outbin']);
	if exist([md.private.runtimename '.tar.gz']) & ~ispc(),
		delete([md.private.runtimename '.tar.gz']);
	end
end

%erase input file if run was carried out on same platform.
hostname=oshostname();
if strcmpi(hostname,cluster.name),
	if md.qmu.isdakota,
		delete([['qmu' num2str(feature('GetPid')) '/'] md.miscellaneous.name '.bin']);
		delete([['qmu' num2str(feature('GetPid')) '/'] md.miscellaneous.name '.queue']);
	else
		delete([md.miscellaneous.name '.bin']);
		delete([md.miscellaneous.name '.toolkits']);
		if ~ispc(),
			delete([md.miscellaneous.name '.queue']);
		else
			delete([md.miscellaneous.name '.bat']);
		end
	end
end
