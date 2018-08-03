function md=parameterize(md,parametername)
%PARAMETERIZE - parameterize a model
%
%   from a parameter matlab file, start filling in all the @model fields that were not 
%   filled in by the mesh.m and mask.m @model methods.
%   Warning: the parameter file must be able to be run in Matlab
%
%   Usage:
%      md=parameterize(md,parametername)
%
%   Example:
%      md=parameterize(md,'Square.par');

%some checks
if ~exist(parametername),
	error(['parameterize error message: file ' parametername ' not found!']);
end

%Try and run parameter file.
temporaryname=['TemporaryParameterFile' num2str(feature('GetPid')) ];
copyfile(parametername,[temporaryname '.m']);

%WARNING: this is a bug of matlab: the TemporaryParameterFile must be cleared
%otherwise matlab keeps the previous version of this file which is not what
%we want!!!
eval(['clear ' temporaryname]);

try,
	eval(temporaryname);
	delete([temporaryname '.m']);
catch me,
	delete([temporaryname '.m']);

	%copy error message
	me2=struct('message',me.message,'stack',me.stack);

	%rename parameter file
	me2.message=regexprep(me2.message,[temporaryname '.m'],parametername);
	for i=1:length(me2.stack)-1,
		me2.stack(i).file=regexprep(me2.stack(i).file,[temporaryname '.m'],parametername);
		me2.stack(i).name=regexprep(me2.stack(i).name,[temporaryname],parametername);
		if strcmp(me2.stack(i).name,'parameterize'),
			%remove error (eval(temporaryname);) misleading
			me2.stack(i)=[];
		end
	end

	%throw error message
	rethrow(me2);
end

%Name and notes
if isempty(md.miscellaneous.name), 
	[path,root,ext]=fileparts(parametername);
	md.miscellaneous.name=root; 
end
md.miscellaneous.notes=['Model created by using parameter file: ' parametername ' on: ' datestr(now)];
