function plot_basaldrag(md,options,width,i,type)

%PLOT_BASALDRAG - plot basal drag
%
%   Usage:
%      plot_basaldrag(md,options,width,i,type);
%
%   See also: PLOTMODEL

%check layer
if dimension(md.mesh)==3,
	if getfieldvalue(options,'layer',1)~=1;
		disp('plot_basaldrag warning: basal drag is displayed in the lower layer')
		changefieldvalue(options,'layer',1);
	end
end

%compute exponents
s=averaging(md,1./md.friction.p,0);
r=averaging(md,md.friction.q./md.friction.p,0);

%compute horizontal velocity
if strcmpi(type,'basal_drag')
	ub=sqrt(md.initialization.vx.^2+md.initialization.vy.^2)/md.constants.yts;
elseif strcmpi(type,'basal_dragx')
	ub=md.initialization.vx/md.constants.yts;
elseif strcmpi(type,'basal_dragy')
	ub=md.initialization.vy/md.constants.yts;
end

%compute basal drag
drag=(max(md.constants.g*(md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base),0)).^r.*(md.friction.coefficient).^2.*ub.^s/1000;

%Figure out if this is a Section plot
if exist(options,'sectionvalue')
	plot_section(md,drag,options,width,i);
	return;
else

	%process data and model
	[x y z elements is2d isplanet]=processmesh(md,[],options);
	[basal_drag datatype]=processdata(md,drag,options);

	%plot basaldrag
	subplot(width,width,i); 
	plot_unit(x,y,z,elements,basal_drag,is2d,isplanet,datatype,options);

	%apply options
	options=addfielddefault(options,'title','Basal drag [kPa]');
	options=addfielddefault(options,'view',2);
	applyoptions(md,basal_drag,options);

end
