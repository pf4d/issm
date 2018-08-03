function [xout,yout] = gdaltransform(x,y,proj_in,proj_out)
%GDALTRANSFORM - switch from one projection system to another 
%
%   Usage:
%      [x,y] = gdaltransform(x1,y1,epsg_in, epsg_out);
%
%   Example: 
%      [x,y] = gdaltransform(md.mesh.long,md.mesh.lat,'EPSG:4326','EPSG:3031')
%
%   For reference: 
%      EPSG: 4326 (lat,long)
%      EPSG: 3411 (Greenland,  UPS 45W, 70N)
%      EPSG: 3031 (Antarctica, UPS 0E,  71S)
%
%      ll2xy default projection Antarctica:
%        +proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs
%      ll2xy default projection Greenland:
%        +proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs
%      Bamber's Greeland projection
%        +proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs
%
%      To get proj.4 string from EPSG, use gdalsrsinfo. Example:    gdalsrsinfo "EPSG:4326" | grep "PROJ.4"

	%give ourselves a unique temporary directory: 
	temproot=tempname; mkdir(temproot);

	fid=fopen([temproot '/.rand1234.txt'],'w');
	fprintf(fid,'%8g %8g\n',[x(:) y(:)]');
	fclose(fid);

	[s,r]=system(['gdaltransform -s_srs ',proj_in,' -t_srs ',proj_out,'  < ' temproot '/.rand1234.txt > ' temproot '/.rand1235.txt']);
	if s~=0 | ~isempty(deblank(r)),
		error(r);
	end
	A=load([temproot '/.rand1235.txt']);
	xout=A(:,1); xout=reshape(xout,size(x));
	yout=A(:,2); yout=reshape(yout,size(y));
