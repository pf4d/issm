function [velx,vely,vel]=sia(md)
%SIA - computation of Shallow Ice velocities
%
%   This routine uses the model of SIA to compute the velocities
%   of a 2d model using the surface slope
%
%   Usage:
%      [velx,vely,vel]=sia(md)

if md.mesh.dimension~=2,
	error('Only 2d meshes are allowed to compute velocity balances');
end

%Get slope
[sx,sy,s]=slope(md);

%Average thickness and B over all elements.
summer=[1;1;1];
hel=md.geometry.thickness(md.mesh.elements)*summer/3;
Bel=md.materials.rheology_B(md.mesh.elements)*summer/3;

Ael=Bel.^(-3);

velx=-2*(md.materials.rho_ice*md.constants.g)^3*s.^2.*sx.*Ael/4.*hel.^4;
vely=-2*(md.materials.rho_ice*md.constants.g)^3*s.^2.*sy.*Ael/4.*hel.^4;
vel=sqrt(velx.^2+vely.^2);
