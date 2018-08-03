function plot_overlay(md,data,options,plotlines,plotcols,i)
%PLOT_OVERLAY - superimpose radar image to a given field
%
%   Usage:
%      plot_overlay(md,data,options,plotlines,plotcols,i)
%
%   See also: PLOTMODEL

%process mesh and data
[x y z elements is2d isplanet]=processmesh(md,[],options);
if strcmpi(data,'none'),
	radaronly=1;
	data=NaN*ones(md.mesh.numberofvertices,1);
	datatype=1;
else
	radaronly=0;
	[data datatype]=processdata(md,data,options);
end

%check is2d
if ~is2d, 
	error('buildoverlay error message: overlay not supported for 3d meshes, project on a layer');
end
if datatype==3,
	error('buildoverlay error message: overlay not supported for quiver plots');
end

%radar power
if ~any(isnan(md.radaroverlay.x)) & ~any(isnan(md.radaroverlay.y)) & ~any(isnan(md.radaroverlay.pwr)),
	disp('plot_overlay info: the radar image held by the model is being used');
	xlim=[min(md.radaroverlay.x) max(md.radaroverlay.x)];
	ylim=[min(md.radaroverlay.y) max(md.radaroverlay.y)];
else
	disp('Extracting radar image...');
	%Get xlim and ylim (used to extract radar image)
	xlim=getfieldvalue(options,'xlim',[min(x) max(x)])/getfieldvalue(options,'unit',1);
	ylim=getfieldvalue(options,'ylim',[min(y) max(y)])/getfieldvalue(options,'unit',1);
	options=addfielddefault(options,'xlim',xlim);
	options=addfielddefault(options,'ylim',ylim);
	md=radarpower(md,options);
end
contrast = getfieldvalue(options,'contrast',1);  
radar = (md.radaroverlay.pwr).^(contrast);
radar = radar./max(radar(:));
if size(radar,3)>1,
	disp('WARNING: color image converted to greyscale intensity image');
	radar=sum(radar,3)/3;
end
%radar(find(radar==0))=1; %Change background from black to white

%InterpFromMeshToGrid
xmin=min(md.radaroverlay.x);
ymax=max(md.radaroverlay.y);
xspacing=(max(md.radaroverlay.x)-min(md.radaroverlay.x))/(length(md.radaroverlay.x));
yspacing=(max(md.radaroverlay.y)-min(md.radaroverlay.y))/(length(md.radaroverlay.y));
nlines=length(md.radaroverlay.y);
ncols =length(md.radaroverlay.x);
disp('Interpolating data on grid...');
if radaronly,
	x_m=xmin:xspacing:xmin+ncols*xspacing;
	y_m=ymax-nlines*yspacing:yspacing:ymax;
	data_grid=NaN*ones(nlines,ncols);
else
	[x_m y_m data_grid]=InterpFromMeshToGrid(elements,x/getfieldvalue(options,'unit',1),y/getfieldvalue(options,'unit',1),...
		data,xmin,ymax,xspacing,yspacing,nlines,ncols,NaN);
end

%Process data_grid (For processing, it is better not to have nan)
pos=find(isinf(data_grid));
if ~isempty(pos),
	disp('Warning: removing Infs from vector (probably log(0)?)');
	data_grid(pos)=NaN;
end
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
data_nan=find(isnan(data_grid));
data_grid(data_nan)=data_min; 

%Special colormaps that require hsv treatment
colorm=getfieldvalue(options,'colormap','Rignot');
if strcmpi(colorm,'Rignot') | strcmpi(colorm,'Seroussi') | strcmpi(colorm,'redblue')
	if strcmpi(colorm,'Rignot'),
		transparency=getfieldvalue(options,'alpha',1);
		h=(data_grid-data_min)/(data_max-data_min+eps);
		if radaronly, h(:)=0; end
		s=max(min((0.1+h).^(1/transparency),1),0);
	elseif strcmpi(colorm,'Seroussi'),
		transparency=getfieldvalue(options,'alpha',1);
		h=1-(data_grid-data_min)/(data_max-data_min+eps)*0.7;
		if radaronly, h(:)=0; end
		s=max(min((0.1+h).^(1/transparency),1),0);
	elseif strcmpi(colorm,'redblue')
		data_mean=data_min+(data_max-data_min)/2;
		h=1*ones(size(data_grid));
		h(find(data_grid<data_mean))=0.7;
		s=max(min(abs(data_grid-data_mean)/(data_max-data_mean) ,1),0);
	else
		error('colormap not supported yet. (''Rignot'' and ''redblue'' are the only cupported colormaps)');
	end
	%(S) Saturation is 0 in NaNs
	s(data_nan)=0;
	%(V) intensity is based on radar image
	v=radar; %use radar power as intensity

	%Transform HSV to RGB
	image_hsv=zeros(size(data_grid,1),size(data_grid,2),3);
	image_hsv(:,:,1)=h; clear h;
	image_hsv(:,:,2)=s; clear s;
	image_hsv(:,:,3)=v; clear v;
	image_rgb=hsv2rgb(image_hsv);
else
	colorm = getcolormap(options);
	len    = size(colorm,1);

	ind = ceil((len-1)*(data_grid-data_min)/(data_max - data_min + eps) +1);
	ind(find(ind>len))=len;
	image_rgb=zeros(size(data_grid,1),size(data_grid,2),3);
	r=colorm(:,1); image_rgb(:,:,1)=r(ind); clear r;
	g=colorm(:,2); image_rgb(:,:,2)=g(ind); clear g;
	b=colorm(:,3); image_rgb(:,:,3)=b(ind); clear b;

	%Now add radarmap
	r = image_rgb(:,:,1).*radar;  r(data_nan) = radar(data_nan);  image_rgb(:,:,1) = r;  clear r;
	g = image_rgb(:,:,2).*radar;  g(data_nan) = radar(data_nan);  image_rgb(:,:,2) = g;  clear g;
	b = image_rgb(:,:,3).*radar;  b(data_nan) = radar(data_nan);  image_rgb(:,:,3) = b;  clear b;
end

%Select plot area 
subplotmodel(plotlines,plotcols,i,options);

%Plot: 
imagesc(md.radaroverlay.x*getfieldvalue(options,'unit',1),md.radaroverlay.y*getfieldvalue(options,'unit',1),image_rgb);set(gca,'YDir','normal');

%last step: mesh overlay?
if exist(options,'edgecolor'),
	hold on
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	patch('Faces',[A B C],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none',...
		'EdgeColor',getfieldvalue(options,'edgecolor'),'LineWidth',getfieldvalue(options,'linewidth',1));
end

%Apply options, without colorbar and without grid
options=changefieldvalue(options,'colormap',colorm);              % We used an HSV colorbar
if ~isnan(data_min),
	options=changefieldvalue(options,'caxis',[data_min data_max]); % force caxis so that the colorbar is ready
end
options=addfielddefault(options,'xlim',xlim);        % default xlim
options=addfielddefault(options,'ylim',ylim);        % default ylim
options=addfielddefault(options,'axis','equal off'); % default axis
applyoptions(md,data,options);
drawnow
