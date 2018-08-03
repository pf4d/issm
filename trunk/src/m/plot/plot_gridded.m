function plot_gridded(md,data,options,plotlines,plotcols,i)
%PLOT_OVERLAY - superimpose radar image to a given field
%
%   Usage:
%      plot_gridded(md,options,plotlines,plotcols,i)
%
%   See also: PLOTMODEL

%process mesh and data
[x y z elements is2d isplanet]=processmesh(md,[],options);
[data datatype]=processdata(md,data,options);

%check is2d
if ~is2d, 
	error('buildgridded error message: gridded not supported for 3d meshes, project on a layer');
end

%Get xlim and ylim (used to extract radar image)
xlim=getfieldvalue(options,'xlim',[min(x) max(x)]);
ylim=getfieldvalue(options,'ylim',[min(y) max(y)]);
postx=getfieldvalue(options,'posting',diff(xlim)/1000);
posty=getfieldvalue(options,'posting',diff(ylim)/1000);

%Interpolating data on grid
[x_m y_m data_grid]=InterpFromMeshToGrid(elements,x,y,data,xlim(1),ylim(2),postx,posty,round(diff(ylim)/posty),round(diff(xlim)/postx),NaN);
data_grid_save = data_grid;
if size(data_grid,1)<3 | size(data_grid,2)<3,
	error('data_grid size too small in plot_gridded, check posting and units');
end

%Process data_grid: add white in NaN and correct caxis accordingly
[data_nani data_nanj]=find(isnan(data_grid) | data_grid==-9999);
if exist(options,'caxis'),
	caxis_opt=getfieldvalue(options,'caxis');
	data_grid(find(data_grid<caxis_opt(1)))=caxis_opt(1);
	data_grid(find(data_grid>caxis_opt(2)))=caxis_opt(2);
	data_min=caxis_opt(1);
	data_max=caxis_opt(2);
else
	data_min=min(data_grid(:));
	data_max=max(data_grid(:));
end

%Select plot area 
subplotmodel(plotlines,plotcols,i,options);

%shading interp;
map    = getcolormap(options);
image_rgb = ind2rgb(uint16((data_grid - data_min)*(length(map)/(data_max-data_min))),map);
if exist(options,'shaded'),
	a    = -45;
	scut = 0.2;
	c    = 1;
	% computes lighting from elevation gradient
	[fx,fy] = gradient(data_grid_save,x_m,y_m);
	fxy = -fx*sind(a) - fy*cosd(a);
	clear fx fy % free some memory...
	fxy(isnan(fxy)) = 0;

	% computes maximum absolute gradient (median-style), normalizes, saturates and duplicates in 3-D matrix
	r = repmat(max(min(fxy/nmedian(abs(fxy),1 - scut/100),1),-1),[1,1,3]);

	% applies contrast using exponent
	rp = (1 - abs(r)).^c;
	image_rgb = image_rgb.*rp;

	% lighter for positive gradient
	k = find(r > 0);
	image_rgb(k) = image_rgb(k) + (1 - rp(k));
end

% set novalues / NaN to black color
if ~isempty(data_nani)
	nancolor=getfieldvalue(options,'nancolor',[1 1 1]);
	image_rgb(sub2ind(size(image_rgb),repmat(data_nani,1,3),repmat(data_nanj,1,3),repmat(1:3,size(data_nani,1),1))) = repmat(nancolor,size(data_nani,1),1);
end

%plot grid
h=imagesc(xlim,ylim,image_rgb);
axis xy

%last step: mesh gridded?
if exist(options,'edgecolor'),
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	patch('Faces',[A B C],'Vertices', [x y z],'FaceVertexCData',data_grid(1)*ones(size(x)),'FaceColor','none','EdgeColor',getfieldvalue(options,'edgecolor'));
end

%Apply options
if ~isnan(data_min) & ~isinf(data_min),
	options=changefieldvalue(options,'caxis',[data_min data_max]); % force caxis so that the colorbar is ready
end
options=addfielddefault(options,'axis','xy equal'); % default axis
applyoptions(md,data,options);

function y = nmedian(x,n)
%NMEDIAN Generalized median filter
%  NMEDIAN(X,N) sorts elemets of X and returns N-th value (N normalized).
%  So:
%     N = 0 is minimum value
%     N = 0.5 is median value
%     N = 1 is maximum value

if nargin < 2
	n = 0.5;
end
y = sort(x(:));
y = interp1(sort(y),n*(length(y)-1) + 1);
