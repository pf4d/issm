function md=ThicknessCorrection(md,varargin)
%THICKNESSCORRECTION - correct the thickness of the ice shelf near the grounding line
%
%   This routine corrects the thickness and the bed on the transition zone
%   by forcing the hydrostatic equilibrium.
%   the thickness is modified as follows:
%      thickness = (1-coeff) * thickness_observation + coeff * thickness_hydrostatic
%   where:
%      coeff=(d/distance)^2;
%      distance=10km by default but can be specified
%
%   Usage:
%      md=ThicknessCorrection(md,varargin);
%
%   Example:
%      md=ThicknessCorrection(md);
%      md=ThicknessCorrection(md,15000);

%initialize thickness with the observations, and get hydrostatic thickness from the dem
thickness=md.geometry.thickness;
thickness_hydro=md.geometry.surface/(1-md.materials.rho_ice/md.materials.rho_water);
hydrostatic_ratio=zeros(size(md.geometry.thickness));

%get nodes on ice sheet and on ice shelf
pos_shelf=find(md.mask.groundedice_levelset<0.);
pos_GL=intersect(unique(md.mesh.elements(find(md.mask.elementongroundedice),:)),unique(md.mesh.elements(find(md.mask.elementonfloatingice),:)));
debug=(length(pos_shelf)>50000);

%check that there is a GL
if isempty(pos_GL)
	error('ThicknessCorrection error message: no grounding line has been detected. Check the model mask');
end

%get distance
if nargin==2,
	distance=varargin{1};
else
	distance=10000;
end

%modify thickness
if (debug), fprintf('%s','      correction progress:   0.00 %'); end
for i=1:length(pos_shelf)

	if (debug & mod(i,100)==0),
		fprintf('\b\b\b\b\b\b\b%5.2f%s',i/length(pos_shelf)*100,' %');
	end

	%search the node on ice sheet the closest to i
	[d posd]=min(sqrt((md.mesh.x(pos_shelf(i))-md.mesh.x(pos_GL)).^2+(md.mesh.y(pos_shelf(i))-md.mesh.y(pos_GL)).^2));

	if d>distance,

		%if d > 15km, hydrostatic equilibrium
		hydrostatic_ratio(pos_shelf(i))=1;
		thickness(pos_shelf(i))=thickness_hydro(pos_shelf(i));

	else

		%else: quadratic combination of hydrostatic equilibrium and observations
		hydrostatic_ratio(pos_shelf(i))=(d/distance)^2;
		thickness(pos_shelf(i))=(1-hydrostatic_ratio(pos_shelf(i)))*thickness(pos_shelf(i))+hydrostatic_ratio(pos_shelf(i))*thickness_hydro(pos_shelf(i));

	end
end
if (debug), fprintf('\b\b\b\b\b\b\b%5.2f%s\n',100,' %'); end

%check the computed thickness
minth=1/(1-md.materials.rho_ice/md.materials.rho_water);
pos=find(isnan(thickness) | (thickness<=0));
thickness(pos)=minth;
hydrostatic_ratio(pos)=-1;

%change bed to take into account the changes in thickness
md.geometry.thickness=thickness;
md.geometry.hydrostatic_ratio=hydrostatic_ratio;
md.geometry.base=md.geometry.surface-md.geometry.thickness;
