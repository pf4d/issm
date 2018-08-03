function expcontourlevelzero(md,mask,level,filename)
%EXPCONTOURLEVELZERO - write an Argus file from a structure recovered from running contourlevelzero 
%
%   Usage:
%      expcontourlevelzero(md,mask,level,filename)
% 
%   Example:
%      expcontourlevelzero(md,md.geometry.thickness,0, 'Level0.exp');
%
%   See also CONTOURLEVELZERO, EXPWRITE

contours=contourlevelzero(md,mask,level);
expwrite(contours,filename);
