function md = transientrestart(md)
%TRANSIENTRESTART - reinitialize model from last transient step
%
%   Usage:
%      md = transientrestart(md)

%Get result and save it again
results = md.results.TransientSolution(end);

newname = ['TransientSolution' num2str(numel(fields(md.results))+1)];
if isfield(md.results,newname)
	error(['Cannot save ' newname ' in md.results']);
else
	disp(['Moving results to ' newname]);
	md.results.(newname) = md.results.TransientSolution;
	md.results.TransientSolution  = struct();
end

%Change time
md.timestepping.start_time = results.time;

%Change initialization fields
if isfield(results,'Vx'),          md.initialization.vx=results.Vx; end
if isfield(results,'Vy'),          md.initialization.vy=results.Vy; end
if isfield(results,'Vz'),          md.initialization.vz=results.Vz; end
if isfield(results,'Temperature'), md.initialization.temperature=results.Temperature; end
if isfield(results,'Pressure'),    md.initialization.pressure=results.Pressure; end

%Deal with new geometry
if isfield(results,'Base') & isfield(results,'Thickness'),
	base=results.Base;
	thickness=results.Thickness;
	if isa(md.mesh,'mesh3dprisms')
		md.mesh.z=base+thickness./md.geometry.thickness.*(md.mesh.z-md.geometry.base);
	end
	md.geometry.base=base;
	md.geometry.thickness=thickness;
	md.geometry.surface=md.geometry.base+md.geometry.thickness;
end

%Update mask
if isfield(results,'MaskGroundediceLevelset'),
	md.mask.groundedice_levelset = results.MaskGroundediceLevelset;
end
if isfield(results,'MaskIceLevelset'),
	md.mask.ice_levelset = results.MaskIceLevelset;
end
