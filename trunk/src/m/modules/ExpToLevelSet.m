function distance = ExpToLevelSet(x,y,contourname);
%EXPTOLEVELSET - Determine levelset distance between a contour and a cloud of points
%
%   Usage:
%      distance=ExpToLevelSet(x,y,contourname);
%
%   x,y:	cloud point
%   contourname:	name of .exp file containing the contours
%   distance:	distance vector representing a levelset where the 0 level is one of the contour segments
%
%   Example:
%      distance=ExpToLevelSet(md.mesh.x,md.mesh.y,'Contour.exp');

% Check usage
if nargin~=3
	help ExpToLevelSet
	error('Wrong usage (see above)');
end

% Call mex module
distance = ExpToLevelSet_matlab(x,y,contourname);
