function [data datatype]=processdata(md,data,options)
%PROCESSDATA - process data to be plotted
%
%   datatype = 1 -> elements
%   datatype = 2 -> nodes
%   datatype = 3 -> node quivers
%   datatype = 4 -> patch
%
%   Usage:
%      [data datatype]=processdata(md,data,options);
%
%   See also: PLOTMODEL, PROCESSMESH

%check format
if (iscell(data) | isempty(data) | length(data)==0 | (length(data)==1 & ~isstruct(data) & isnan(data))),
	error('plotmodel error message: data provided is empty');
end

%specials for struct
if isstruct(data),
	disp('data provided is a struct with the following fields:');
	F=fields(data);
	for i=1:numel(F),
		disp(['   ' num2str(i) ': ' F{i} ]);
	end
	choice=input(['please enter the field number? (between 1 and ' num2str(numel(F)) ')  ']);
	[data datatype]=processdata(md,data(1).(F{choice}),options);
end

%Process NaN if any (do not know before mask is applied)
if exist(options,'nan')
	data(find(isnan(data)))=getfieldvalue(options,'nan',0);
end

%special case for mesh 2dvertical
if strcmp(domaintype(md.mesh),'2Dvertical'),
	[data datatype] = processdata(md.mesh,md,data,options);
	return;
end

%needed later on
if isprop(md.mesh,'numberofvertices2d'), 
	numberofvertices2d=md.mesh.numberofvertices2d; 
	numberofelements2d=md.mesh.numberofelements2d; 
else 
	numberofvertices2d=NaN;
	numberofelements2d=NaN;
end

%initialize datatype
datatype=0;

%get datasize
datasize=size(data);

%transpose data if necessary
if (size(data,2) > size(data,1)),
	data=data';
end
datasize=size(data);

%convert to double if necessary
if ~isnumeric(data);
	disp('processdata info message: data is not numeric (logical?). Converted to double');
	data=double(data);
end

%check length
if datasize(1)~=md.mesh.numberofvertices & datasize(1)~=md.mesh.numberofelements & datasize(1)~=md.mesh.numberofvertices*6 & (strcmp(md.mesh.domaintype(),'3D') & ~(datasize(1)==numberofelements2d | datasize(1)==numberofvertices2d))
	error('plotmodel error message: data not supported yet');
end

%quiver?
if datasize(2)>1,
	datatype=3;

	%check number of columns, add zeros if necessary,
	if (dimension(md.mesh)==3)
		if datasize(2)==2,
			data=[data, zeros(datasize(1),1)];
		elseif datasize(2)~=3,
			error('plotmodel error message: data provided should have 2 or 3 columns for quiver plot, and 1 for regular plot');
		end
	end
end

%treat the case datasize(1)=6*nodes
if datasize(1)==6*md.mesh.numberofvertices
	%keep the only norm of data
	data1=data(1:6:md.mesh.numberofvertices*6,:);
	data2=data(2:6:md.mesh.numberofvertices*6,:);
	data=sqrt(data1.^2+data2.^2);
	datasize(1)=md.mesh.numberofvertices;
	%---> go to node data
end

%treat the case datasize(1)=nodes2d
if (dimension(md.mesh)==3 & datasize(1)==numberofvertices2d),
	data=project3d(md,'vector',data,'type','node');
	datasize(1)=md.mesh.numberofvertices;
	%---> go to node data
end

%treat the case datasize(1)=nodes2d
if (dimension(md.mesh)==3 & datasize(1)==numberofelements2d),
	data=project3d(md,'vector',data,'type','element');
	datasize(1)=md.mesh.numberofelements;
	%---> go to node data
end

%smoothing?
if exist(options,'smooth')
	data=averaging(md,data,getfieldvalue(options,'smooth'));
	datasize(1)=md.mesh.numberofvertices;
	%---> go to node data
end

%element data
if (datasize(1)==md.mesh.numberofelements & datasize(2)==1),

	%Initialize datatype if non patch
	if datatype~=4 & datatype~=5,
		datatype=1;
	end

	%Mask?
	if exist(options,'mask'),
		flags=getfieldvalue(options,'mask');
		maskvalue=getfieldvalue(options,'maskvalue',NaN);
		pos=find(~flags);
		if length(flags)==md.mesh.numberofvertices,
			[pos2 dummy]=find(ismember(md.mesh.elements,pos));
			data(pos2,:)=maskvalue;
		elseif length(flags)==md.mesh.numberofelements
			data(pos,:)=maskvalue;
		else
			disp('plotmodel warning: mask length not supported yet (supported length are md.mesh.numberofvertices and md.mesh.numberofelements');
		end
	end

	%log?
	if exist(options,'log'),
		bounds=getfieldvalue(options,'caxis',[min(data(:)) max(data(:))]);
		data(find(data<bounds(1)))=bounds(1);
		if any(data<=0),
			error('Log option cannot be applied on negative values. Use caxis option (Rignot''s settings: [1.5 max(data)])');
		end
		pos=find(~isnan(data));
		data(pos)=log(data(pos))/log(getfieldvalue(options,'log'));
	end
end

%node data
if (datasize(1)==md.mesh.numberofvertices & datasize(2)==1),
	datatype=2;

	%Mask?
	if exist(options,'mask'),
		flags=getfieldvalue(options,'mask');
		maskvalue=getfieldvalue(options,'maskvalue',NaN);
		pos=find(~flags);
		if length(flags)==md.mesh.numberofvertices,
			data(pos,:)=maskvalue;
		elseif length(flags)==md.mesh.numberofelements
			data(md.mesh.elements(pos,:),:)=maskvalue;
		else
			disp('plotmodel warning: mask length not supported yet (supported length are md.mesh.numberofvertices and md.mesh.numberofelements');
		end
	end

	%log?
	if exist(options,'log'),
		bounds=getfieldvalue(options,'caxis_pre',[min(data(:)) max(data(:))]);
		data(find(data<bounds(1)))=bounds(1);
		if any(data<=0),
			error('Log option cannot be applied on negative values. Use caxis option (Rignot''s settings: [1.5 max(data)])');
		end
		pos=find(~isnan(data));
		data(pos)=log(data(pos))/log(getfieldvalue(options,'log'));
	end
end

%layer projection? 
if getfieldvalue(options,'layer',0)>=1,
	data=project2d(md,data,getfieldvalue(options,'layer')); %project onto 2d mesh
end

%control arrow density if quiverplot
if datatype==3 & exist(options,'density')
	databak=data;
	data=NaN*ones(datasize);
	density=getfieldvalue(options,'density');
	data(1:density:end,:)=databak(1:density:end,:);
	clear databak
end
if datatype==3,
	%Mask?
	if exist(options,'mask'),
		flags=getfieldvalue(options,'mask');
		pos=find(~flags);
		if length(flags)==md.mesh.numberofvertices,
			data(pos,:)=NaN;
		elseif length(flags)==md.mesh.numberofelements
			data(md.mesh.elements(pos,:),:)=NaN;
		else
			disp('plotmodel warning: mask length not supported yet (supported length are md.mesh.numberofvertices and md.mesh.numberofelements');
		end
	end
end

%OK, if datatype=0 error out
if datatype==0,
	error(['data provided not recognized or not supported']);
end
