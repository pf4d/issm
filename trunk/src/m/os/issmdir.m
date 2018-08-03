function ISSM_DIR=issmdir()
%ISSMDIR - Get ISSM_DIR environment variable
%
%   Usage:
%      ISSM_DIR=issmdir()

%Initialize output ISSM_DIR
ISSM_DIR='';
slash=filesep();

%Get ISSM_DIR from function path (we do not want to force users to edit their bashrc)
path=which('issmdir');

%issmdir might be in bin,
pos=strfind(path,['bin' slash 'issmdir.m']);
if ~isempty(pos),
	ISSM_DIR=path(1:pos-1);
else
	pos=strfind(path,['src' slash 'm' slash 'os' slash 'issmdir.m']);
	if ~isempty(pos),
		ISSM_DIR=path(1:pos-1);
	end
end

if isempty(ISSM_DIR),
	error('Could not determine the location of ISSM...');
end
