function mask = gmtmask(lat,long,varargin)
%GMTMASK - figure out which lat,long points are on the ocean
%
%   Usage:
%      mask.ocean = gmtmask(md.mesh.lat,md.mesh.long);
%

	%are we doing a recursive call? 
	if nargin==3,
		recursive=1;
	else 
		recursive=0;
	end
	
	if(recursive)disp(sprintf('             recursing: num vertices %i',length(lat)));
	else disp(sprintf('gmtmask: num vertices %i',length(lat)));
	end
	
	%Check lat and long size is not more than 50,000; If so, recursively call gmtmask: 
	if length(lat)>50000,
		for i=1:50000:length(lat),
			j=i+50000-1;
			if j>length(lat),
				j=length(lat);
			end
			mask(i:j)=gmtmask(lat(i:j),long(i:j),1);
		end
		return
	end
	
	%First, write our lat,long file for gmt:
	nv=length(lat);
	dlmwrite('./all_vertices.txt',[long lat (1:nv)'],'delimiter','\t','precision',10);

	%Avoid bypassing of the ld library path by Matlab (:()
	if ismac,
		dyld_library_path_old=getenv('DYLD_LIBRARY_PATH');
		setenv('DYLD_LIBRARY_PATH',[ issmdir '/externalpackages/curl/install/lib:' issmdir '/externalpackages/hdf5/install/lib:' issmdir '/externalpackages/netcdf/install/lib' ]);
	end

	%figure out which vertices are on the ocean, which one on the continent:
	[status,result] = system('gmt gmtselect ./all_vertices.txt -h0 -Df -R0/360/-90/90  -A0 -JQ180/200 -Nk/s/s/k/s > ./oce_vertices.txt');
	if status~=0,
		error(result);
	end

	%reset DYLD_LIBRARY_PATH to what it was: 
	if ismac,
		setenv('DYLD_LIBRARY_PATH',dyld_library_path_old);
	end
	%read the con_vertices.txt file and flag our mesh vertices on the continent
	fid=fopen('./oce_vertices.txt','r');
	line=fgets(fid); 
	line=fgets(fid);
	oce_vertices=[];
	while line~=-1,
		ind=str2num(line); ind=ind(3);
		oce_vertices=[oce_vertices;ind];
		line=fgets(fid);
	end


	mask=zeros(nv,1);
	mask(oce_vertices)=1;
	
	system('rm -rf ./all_vertices.txt ./oce_vertices.txt ./gmt.history');
	if ~recursive, disp(sprintf('gmtmask: done')); end;
