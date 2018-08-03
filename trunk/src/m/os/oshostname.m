function hostname=oshostname()
%OSHOSTNAME - Determine hostname, irrespective of os type
%
%   Usage:
%      hostname=oshostname();

%See http://www.mathworks.com/help/matlab/ref/system.html "tips" section
%We need to add < /dev/null otherwise what is in the clipboard is added
[status,hostname]=system('hostname < /dev/null');

% If that command did not work, we have an alternative
if status~=0,
	if ispc
		hostname = getenv('COMPUTERNAME');
	else
		hostname = getenv('HOSTNAME');
	end
end

% Take out minus signs
hostname = strrep(hostname,'-','');

% Trim and lower case
hostname = strtrim(lower(hostname));

% Check that machine name is not empty
if isempty(hostname),
	error('Cannot determine machine name');
end
