function md=sethydrostaticmask(md)
%SETHYDROSTATICMASK - establish groundedice_levelset field
%
%   Determines grounded and floating ice position based on 
%   md.geometry.bed and md.geometry.thickness
%
%   Usage:
%      md=sethydrostaticmask(md)
%
%   Examples:
%      md=sethydrostaticmask(md);

if(length(md.geometry.bed)~=md.mesh.numberofvertices | length(md.geometry.thickness)~=md.mesh.numberofvertices | length(md.geometry.base)~=md.mesh.numberofvertices),
		error('hydrostaticmask error message: fields in md.geometry do not have the right size.');
end

%grounded ice level set
md.mask.groundedice_levelset=md.geometry.thickness+md.geometry.bed*md.materials.rho_water/md.materials.rho_ice;

%Check consistency of geometry
pos=find(md.mask.groundedice_levelset>0);
if(any(md.geometry.base(pos)~=md.geometry.bed(pos))),
	disp('WARNING: md.geometry.bed and md.geometry.base not equal on grounded ice');
end

pos=find(md.mask.groundedice_levelset<=0);
if(any(md.geometry.base(pos)<md.geometry.bed(pos))),
	disp('WARNING: md.geometry.base < md.geometry.bed on floating ice');
end
