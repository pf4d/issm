function [bx by b]=basalstress(md)
%BASALSTRESS - compute basal stress from basal drag and geometric information. 
%
%      Computes basal stress from geometric information and ice velocity in md.initialization.
%
%   Usage:
%      [bx by b]=basalstress(md);
%
%   See also: plot_basaldrag

%compute exponents
s=averaging(md,1./md.friction.p,0);
r=averaging(md,md.friction.q./md.friction.p,0);

%compute horizontal velocity
ub=sqrt(md.initialization.vx.^2+md.initialization.vy.^2)/md.constants.yts;
ubx=md.initialization.vx/md.constants.yts;
uby=md.initialization.vy/md.constants.yts;

%compute basal drag
bx=(md.constants.g*(md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base)).^r.*(md.friction.coefficient).^2.*ubx.^s;
by=(md.constants.g*(md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base)).^r.*(md.friction.coefficient).^2.*uby.^s;
b=(md.constants.g*(md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base)).^r.*(md.friction.coefficient).^2.*ub.^s;
