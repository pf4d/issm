function Struct=shpread(filename)
%SHPREAD - read a shape file and build a Structure
%
%   This routine reads a shape file .shp and builds a Structure containing the 
%   fields x and y corresponding to the coordinates, one for the filename of
%   the shp file, for the density, for the nodes, and a field closed to 
%   indicate if the domain is closed. 
%   If this initial shapefile is point only, the fields closed and
%   points are ommited
%   The first argument is the .shp file to be read and the second one (optional) 
%   indicates if the last point shall be read (1 to read it, 0 not to).
%
%   Usage:
%      Struct=shpread(filename)
%
%   Example:
%      Struct=shpread('domainoutline.shp')
%
%   See also EXPDOC, EXPWRITEASVERTICES

%some checks
if ~exist(filename),
	error(['shpread error message: file ' filename ' not found!']);
end

%initialize number of profile
count=0;

%read shapefile
shp=shaperead(filename);

Struct=struct([]);
fields=fieldnames(shp);
for i=1:length(shp),
	if strcmpi(shp(i).Geometry,'Polygon'),
		x=shp(i).X'; y=shp(i).Y';
		ids=find(isnan(x));
		x(ids)=[]; y(ids)=[];

		Struct(end+1).x=x;
		Struct(end).y=y;
		Struct(end).nods=length(x);
		Struct(end).density=1;
		Struct(end).closed=1;
		if isfield(shp,'id'),
			Struct(end).name=num2str(shp(i).id);
		else
			Struct(end).name='';
		end
		for j=1:length(fields),
			field=fields{j};
			if ~(strcmpi(field,'X') | strcmpi(field,'Y') | strcmpi(field,'id')),
				Struct(end).(field)=shp(i).(field);
			end
		end
	end
	
	if strcmpi(shp(i).Geometry,'Line'),
		x=shp(i).X'; y=shp(i).Y';
		ids=find(isnan(x));
		x(ids)=[]; y(ids)=[];

		Struct(end+1).x=x;
		Struct(end).y=y;
		Struct(end).nods=length(x);
		Struct(end).density=1;
		Struct(end).closed=1;
		if isfield(shp,'id'),
			Struct(end).name=num2str(shp(i).id);
		else
			Struct(end).name='';
		end
		for j=1:length(fields),
			field=fields{j};
			if ~(strcmpi(field,'X') | strcmpi(field,'Y') | strcmpi(field,'id')),
				Struct(end).(field)=shp(i).(field);
			end
		end
	end


	if strcmpi(shp(i).Geometry,'Point'),
		x=shp(i).X'; y=shp(i).Y';
		ids=find(isnan(x));
		x(ids)=[]; y(ids)=[];

		Struct(end+1).x=x;
		Struct(end).y=y;
		Struct(end).density=1;
		if isfield(shp,'id'),
			Struct(end).name=num2str(shp(i).id);
		else
			Struct(end).name='';
		end
		for j=1:length(fields),
			field=fields{j};
			if ~(strcmpi(field,'X') | strcmpi(field,'Y') | strcmpi(field,'id')),
				Struct(end).(field)=shp(i).(field);
			end
		end
	end
end

end
