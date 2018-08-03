function qmumarshall(md.qmu.variables,responses)
%QMUMARSHALL - output ISSM compatible binary file with qmu fields. This is 
%   in addition to the marshall routine for regular solve routines.
%   Usage:
%      qmumarshall(md.qmu.variables,responses)
% 
%   where variables and responses are the Dakota variables and responses found in the model @md.

%some checks on list of arguments
if ((nargin~=3) & (nargout~=0))
	qmumarshallusage;
	error('marshall error message');
end

disp(['qmu marshalling file ' md.miscellaneous.name '.bin']);

%open file for binary adding 
fid=fopen([ md.miscellaneous.name '.bin'],'ab');
if fid==-1,
	error(['qmumarshall error message: could not open ' [md.miscellaneous.name '.bin'],' file for binary adding']);
end

%deal with variables
WriteData(fid,md.numvariabledescriptors,'Integer','numvariabledescriptors');
for i=1:md.numvariabledescriptors,
	field_name=md.qmu.variabledescriptors{i};
	WriteData(fid,field_name,'String',['variabledescriptor' num2str(i)]);
end

%deal with responses

%write number of responses to disk
WriteData(fid,md.qmu.numberofresponses,'Integer','numberofresponses');
WriteData(fid,md.qmu.numresponsedescriptors,'Integer','numresponsedescriptors');
for i=1:md.qmu.numresponsedescriptors,
	field_name=md.qmu.responsedescriptors{i};
	WriteData(fid,field_name,'String',['responsedescriptor' num2str(i)]);
end

%write response specific data
qmu_segments=0;

for i=1:numel(md.qmu.responsedescriptors),
	field_name=md.qmu.responsedescriptors{i};
	if strncmpi(field_name,'indexed_MassFlux',16),
		qmu_segments=1;
	end
end

if qmu_segments,
	WriteData(fid,md.qmu.mass_flux_num_profiles,'Integer','qmu_mass_flux_num_profiles');
	for i=1:md.qmu.mass_flux_num_profiles,
		WriteData(fid,md.qmu.mass_flux_segments{i},'Mat',['qmu_mass_flux_segments' num2str(i)]);
	end
else
	md.qmu.mass_flux_num_profiles=0;
	WriteData(fid,md.qmu.mass_flux_num_profiles,'Integer','qmu_mass_flux_num_profiles');
end

%write part and npart to disk
WriteData(fid,md.qmu.numberofpartitions,'Integer','npart');
WriteData(fid,md.qmu.partition,'Mat','part');

%close file
st=fclose(fid);
if st==-1,
	error(['qmumarshall error message: could not close file ' [md.miscellaneous.name '.bin']]);
end

end

function qmumarshallusage()
disp(' ');
disp('function qmumarshall(md.qmu.variables,responses)');
end
