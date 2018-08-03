function md=preqmu(md,options)
%QMU - apply Quantification of Margins and Uncertainties techniques 
%      to a solution sequence (like stressbalance.m, progonstic.m, etc ...), 
%      using the Dakota software from Sandia.
%
%   options come from the solve.m routine. They can include Dakota options:
%
%       qmudir:  any directory where to run the qmu analysis
%       qmufile: input file for Dakota
%       ivar: selection number for variables input (if several are specified in variables)
%       iresp: same thing for response functions
%       imethod: same thing for methods
%       iparams: same thing for params
%       overwrite: overwrite qmudir before analysis

disp('preprocessing dakota inputs');
qmudir    = getfieldvalue(options,'qmudir',['qmu' num2str(feature('GetPid'))]);  % qmudir = ['qmu_' datestr(now,'yyyymmdd_HHMMSS')];
qmufile   = getfieldvalue(options,'qmufile','qmu');% qmufile cannot be changed unless ????script.sh is also changed
overwrite = getfieldvalue(options,'overwrite','n');
ivar      = getfieldvalue(options,'ivar',1);
iresp     = getfieldvalue(options,'iresp',1);
imethod   = getfieldvalue(options,'imethod',1);
iparams   = getfieldvalue(options,'iparams',1);

%first create temporary directory in which we will work
if strncmpi(overwrite,'y',1)
	system(['rm -rf ' qmudir '/*']); 
else
	%does the directory exist? if so, then error out
	if exist(qmudir)==7,
		error('Existing ''%s'' directory, cannot overwrite. Specify ''overwrite'',''y'' option in solve arguments.',options.qmudir);
	end
end
mkdir(qmudir)
cd(qmudir)

%when running in library mode, the in file needs to be called md.miscellaneous.name.qmu.in
qmufile=[md.miscellaneous.name ];

%retrieve variables and resposnes for this particular analysis.
variables=md.qmu.variables(ivar);
responses=md.qmu.responses(iresp);

%expand variables and responses
variables=expandvariables(md,variables);
responses=expandresponses(md,responses);

%go through variables and responses, and check they don't have more than md.qmu.numberofpartitions values. Also determine numvariables and numresponses{{{
numvariables=0;
variable_fieldnames=fieldnames(variables);
for i=1:length(variable_fieldnames),
	field_name=variable_fieldnames{i};
	fieldvariables=variables.(field_name);
	for j=1:numel(fieldvariables)
		if strncmpi(fieldvariables(j).descriptor,'scaled_',7) && str2int(fieldvariables(j).descriptor,'last')>md.qmu.numberofpartitions,
			error('preqmu error message: one of the expanded variables has more values than the number of partitions (setup in md.qmu.numberofpartitions)');
		end
	end
	numvariables=numvariables+numel(variables.(field_name));
end

numresponses=0;
response_fieldnames=fieldnames(responses);
for i=1:length(response_fieldnames),
	field_name=response_fieldnames{i};
	fieldresponses=responses.(field_name);
	for j=1:numel(fieldresponses)
		if strncmpi(fieldresponses(j).descriptor,'scaled_',7) && str2int(fieldresponses(j).descriptor,'last')>md.qmu.numberofpartitions,
			error('preqmu error message: one of the expanded responses has more values than the number of partitions (setup in md.qmu.numberofpartitions)');
		end
	end
	numresponses=numresponses+numel(responses.(field_name));
end
%}}}}

%create in file for dakota
dakota_in_data(md.qmu.method(imethod),variables,responses,md.qmu.params(iparams),qmufile);
system(['rm -rf ' md.miscellaneous.name '.m']);

%build a list of variables and responses descriptors. the list is not expanded. {{{
variabledescriptors={};
variable_fieldnames=fieldnames(md.qmu.variables(ivar));
for i=1:length(variable_fieldnames),
	field_name=variable_fieldnames{i};
	fieldvariables=md.qmu.variables(ivar).(field_name);
	for j=1:numel(fieldvariables)
		variabledescriptors{end+1}=fieldvariables(j).descriptor;
	end
end

responsedescriptors={};
response_fieldnames=fieldnames(md.qmu.responses(iresp));
for i=1:length(response_fieldnames),
	field_name=response_fieldnames{i};
	fieldresponses=md.qmu.responses(iresp).(field_name);
	for j=1:numel(fieldresponses)
		responsedescriptors{end+1}=fieldresponses(j).descriptor;
	end
end
%}}}

%register the fields that will be needed by the Qmu model.
md.qmu.numberofresponses=numresponses;
md.qmu.variabledescriptors=variabledescriptors;
md.qmu.responsedescriptors=responsedescriptors;

%now, we have to provide all the info necessary for the solutions to compute the responses. For ex, if mass_flux 
%is a response, we need a profile of points.  For a misfit, we need the observed velocity, etc ...
md=process_qmu_response_data(md);
