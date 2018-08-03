function outoptions=process_qmu_options(options)
%PROCESS_QMU_OPTIONS - set up default options for qmu phase
%
%   Usage:
%      options=process_qmu_options(options)
%
%   See also: QMU,RECOVER_QMU_OPTIONS

%analysis_type: check on this option, error out otherwise
found=0;
for i=1:size(options,1),
	if strcmpi(options{i,1},'analysis_type'),
		analysis_type=options{i,2};
		found=1;
	end
end
if ~found,
	error('recover_qmu_options error message: no ''analysis_type'' was provided');
end

%package: is there one? default to ''JPL''
found=0;
for i=1:size(options,1),
	if strcmpi(options{i,1},'package'),
		package=options{i,2};
		found=1;
	end
end
if ~found,
	disp('recover_qmu_options info message: no ''package'' was provided, defaulting to ''JPL''');
	options(end+1,:)={'package' 'JPL'};
	package='JPL';
end

if ~ischar(package), 
	error(['process_qmu_options error message: package ' package ' not supported yet']);
end

%check solution type is supported
if ~(strcmpi(analysis_type,'control') |  ...
		strcmpi(analysis_type,'stressbalance') |  ...
		strcmpi(analysis_type,'masstransport') |  ...
		strcmpi(analysis_type,'thermal') |  ...
		strcmpi(analysis_type,'parameters') |  ...
		strcmpi(analysis_type,'transient') ),
	error(['process_qmu_options error message: analysis_type ' analysis_type ' not supported yet!']);
end

%  process qmu arguments

%first, the defaults
qmudir ='qmu';% qmudir =['qmu_' datestr(now,'yyyymmdd_HHMMSS')];
qmufile='qmu';
ivar   =1;
iresp  =1;
imethod=1;
iparams=1;
runmpi =false;

for i=1:size(options,1),
	switch options{i,1},
	case 'qmudir'
		qmudir=options{i,2};
	case 'qmufile'
		qmufile=options{i,2};
	case 'ivar'
		ivar=options{i,2};
	case 'iresp'
		iresp=options{i,2};
	case 'imethod'
		imethod=options{i,2};
	case 'iparams'
		iparams=options{i,2};
	case 'overwrite'
		outoptions.overwrite=options{i,2};
	case 'keep'
		outoptions.keep=options{i,2};
	case 'outfiles'
		outoptions.outfiles=options{i,2};
	case 'rstfile'
		outoptions.rstfile=options{i,2}; 
	case 'rundakota'
		outoptions.rundakota=options{i,2};
	case 'runmpi'
		runmpi=options{i,2};
	otherwise
		%nothing
	end
end

%setup final options structure
outoptions.analysis_type=analysis_type;
outoptions.package=package;
outoptions.qmudir=qmudir;
outoptions.qmufile=qmufile;
outoptions.ivar=ivar;
outoptions.iresp=iresp;
outoptions.imethod=imethod;
outoptions.iparams=iparams;
outoptions.runmpi=runmpi;
