%STARTUP - Matlab startup script
%
%   startup.m is a script run by matlab at the beginning of a session, just
%   before handing over the prompt to the user. This delivery startup.m script
%   has been customized here for the ISSM code. This startup script should be
%   run by users before trying to use ISSM. The best way to do that is to put
%   the startup file in the location where Matlab starts and established its
%   root directory.

% clear the last warning to focus on the warnings of the ISSM path
lastwarn(''); 

%Recover ISSM_TIER , or if on a Windows machine, ISSM_TIER_WIN
ISSM_TIER=getenv('ISSM_TIER_WIN');

if (isempty(ISSM_TIER)),
	error('issmdir error message: ''ISSM_TIER'' environment variable is empty! You should define ISSM_TIER in your .cshrc or .bashrc!');
end

%Now add all issm code paths necessary to run issm smoothly. 
%We capture the error output, so that we can warn the user to update 
%the variable ISSM_TIER in this file, in case it is not correctly setup. 

%ISSM path
addpath(pwd); %add current path first
addpath([pwd '\bin']);
addpath([pwd '\lib']);

%Check on any warning messages that might indicate that the paths were not correct. 
if ~isempty(lastwarn),
	fprintf('\n  Error trying to setup ''ISSM'' code paths. Try and update the ISSM_TIER variable in your .cshrc or .bashrc!\n');
	fprintf('  ''ISSM'' will not  work at all until this is resolved\n\n');
else
	fprintf('\n  To get started with ISSM, type issmdoc at the command prompt.\n\n');
end

%disable matlab bell!
beep off;

% no warning if we try to plot while in nojvm (will not be supported in future releases)
warning off MATLAB:HandleGraphics:noJVM

%at the end, get to tests directory if ISSM_TESTS exists: 
ISSM_TESTS=getenv('ISSM_TESTS');
if ~isempty(ISSM_TESTS),
	cd(ISSM_TESTS);
end
