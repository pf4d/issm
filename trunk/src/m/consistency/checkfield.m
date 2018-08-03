function md = checkfield(md,varargin)
%CHECKFIELD - check field consistency
%
%   Used to check model consistency.
%   Requires: 
%     'field' or 'fieldname' option. If 'fieldname' is provided, it will retrieve it from the model md. (md.(fieldname)) 
%             If 'field' is provided, it will assume the argument following 'field' is a numeric array.
%   Available options:
%      - NaN: 1 if check that there is no NaN
%      - Inf: 1 if check that there is no Inf
%      - size: [lines cols], NaN for non checked dimensions
%      - >:  greater than provided value
%      - >=: greater or equal to provided value
%      - <:  smallerthan provided value
%      - <=: smaller or equal to provided value
%      - < vec:  smallerthan provided values on each vertex
%      - timeseries: 1 if check time series consistency (size and time)
%      - values: cell of strings or vector of acceptable values
%      - numel: list of acceptable number of elements
%      - cell: 1 if check that is cell
%      - empty: 1 if check that non empty
%      - message: overloaded error message
%
%   Usage:
%      md = checkfield(md,fieldname,options);

%get options
options=pairoptions(varargin{:});

%get field: 
if exist(options,'field'), 
	field=getfieldvalue(options,'field'); 
	fieldname=getfieldvalue(options,'fieldname','no fieldname'); 
else
	fieldname=getfieldvalue(options,'fieldname'); 
	eval(['field=md.' fieldname ';']);
end

%check empty
if exist(options,'empty')
	if isempty(field),
		md = checkmessage(md,getfieldvalue(options,'message',...
			['field ''' fieldname ''' is empty']));
	end
end

%Check size
if exist(options,'size')
	fieldsize=getfieldvalue(options,'size');
	if isnan(fieldsize(1)),
		if (size(field,2)~=fieldsize(2)),
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' should have ' num2str(fieldsize(2)) ' columns']));
		end
	elseif isnan(fieldsize(2)),
		if (size(field,1)~=fieldsize(1)),
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' should have ' num2str(fieldsize(1)) ' lines']));
		end
	else
		if ((size(field,1)~=fieldsize(1)) |  (size(field,2)~=fieldsize(2)))
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' size should be ' num2str(fieldsize(1)) ' x ' num2str(fieldsize(2))]));
		end
	end
end

%Check numel
if exist(options,'numel')
	fieldnumel=getfieldvalue(options,'numel');
	if ~ismember(numel(field),fieldnumel),
		if length(fieldnumel)==1
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' size should be ' sprintf('%g ',fieldnumel) ]));
		elseif length(fieldnumel)==2
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' size should be ' num2str(fieldnumel(1)) ' or ' num2str(fieldnumel(2)) ]));
		else
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' size should be ' sprintf('%g, ',fieldnumel(1:end-1)) ' or ' num2str(fieldnumel(end)) ]));
		end
	end
end

%check NaN
if getfieldvalue(options,'NaN',0);
	field2=reshape(field,prod(size(field)),1);
	if any(isnan(field2)),
		md = checkmessage(md,getfieldvalue(options,'message',...
			['NaN values found in field ''' fieldname '''']));
	end
end

%check Inf
if getfieldvalue(options,'Inf',0);
	field2=reshape(field,prod(size(field)),1);
	if any(isinf(field2)),
		md = checkmessage(md,getfieldvalue(options,'message',...
			['Inf values found in field ''' fieldname '''']));
	end
end


%check cell
if getfieldvalue(options,'cell',0);
	if ~iscell(field),
		md = checkmessage(md,getfieldvalue(options,'message',...
			['field ''' fieldname ''' should be a cell']));
	end
end

%check values
if exist(options,'values')
	fieldvalues=getfieldvalue(options,'values');
	if iscell(fieldvalues), %strings
		if ischar(field) | iscell(fieldvalues),
			if any(~ismember(field,fieldvalues)),
				if length(fieldvalues)==1
					md = checkmessage(md,getfieldvalue(options,'message',...
						['field ''' fieldname ''' value should be ''' fieldvalues{1} '''']));
				elseif length(fieldvalues)==2
					md = checkmessage(md,getfieldvalue(options,'message',...
						['field ''' fieldname ''' values should be ''' fieldvalues{1} ''' or ''' fieldvalues{2} '''']));
				else
					md = checkmessage(md,getfieldvalue(options,'message',...
						['field ''' fieldname ''' should have values in ' sprintf('''%s'', ',fieldvalues{1:end-1}) 'or ''' fieldvalues{end} '''']));
				end
			end
		else
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' should be one of the following strings: ' sprintf('''%s'', ',fieldvalues{1:end-1}) 'or ''' fieldvalues{end} '''']));
		end
	else
		field2=reshape(field,prod(size(field)),1);
		if isnumeric(field),
			if any(~ismember(field2,fieldvalues)),
				md = checkmessage(md,getfieldvalue(options,'message',...
					['field ''' fieldname ''' should have values in [' num2str(fieldvalues) ']']));
			end
		else
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' should be a number in [' num2str(fieldvalues) ']']));
		end
	end
end

%check greater
if exist(options,'>=')
	lowerbound=getfieldvalue(options,'>=');
	field2=reshape(field,prod(size(field)),1);
	if getfieldvalue(options,'timeseries',0), field2=reshape(field(1:end-1,:),prod(size(field(1:end-1,:))),1); end
	if any(field2<lowerbound),
		md = checkmessage(md,getfieldvalue(options,'message',...
			['field ''' fieldname ''' should have values above ' num2str(lowerbound(1,1))]));
	end
end
if exist(options,'>')
	lowerbound=getfieldvalue(options,'>');
	field2=reshape(field,prod(size(field)),1);
	if getfieldvalue(options,'timeseries',0), field2=reshape(field(1:end-1,:),prod(size(field(1:end-1,:))),1); end
	if any(field2<=lowerbound),
		md = checkmessage(md,getfieldvalue(options,'message',...
			['field ''' fieldname ''' should have values above ' num2str(lowerbound(1,1))]));
	end
end

%check smaller
if exist(options,'<=')
	upperbound=getfieldvalue(options,'<=');
	field2=reshape(field,prod(size(field)),1);
	if getfieldvalue(options,'timeseries',0), field2=reshape(field(1:end-1,:),prod(size(field(1:end-1,:))),1); end
	if any(field2>upperbound),
		md = checkmessage(md,getfieldvalue(options,'message',...
			['field ''' fieldname ''' should have values below ' num2str(upperbound(1,1))]));
	end
end
if exist(options,'<')
	upperbound=getfieldvalue(options,'<');
	field2=reshape(field,prod(size(field)),1);
	if getfieldvalue(options,'timeseries',0), field2=reshape(field(1:end-1,:),prod(size(field(1:end-1,:))),1); end
	if any(field2>=upperbound),
		md = checkmessage(md,getfieldvalue(options,'message',...
			['field ''' fieldname ''' should have values below ' num2str(upperbound(1,1))]));
	end
end

%Check row of stringrow
if getfieldvalue(options,'stringrow',0),
	if(size(field,1)~=1 & size(field,1)~=0),
		md = checkmessage(md,getfieldvalue(options,'message',...
			['field ''' fieldname ''' should have only one row']));
	end
	if ~iscell(field),
		md = checkmessage(md,getfieldvalue(options,'message',...
			['field ''' fieldname ''' should be a cell of strings']));
	else
		for i=1:size(field,2),
			if ~ischar(field{i}),
				md = checkmessage(md,getfieldvalue(options,'message',...
					['field ''' fieldname ''' values should a cell of chars']));
			end
		end
	end
end

%check file
if getfieldvalue(options,'file',0),
	if ~exist(field,'file')
		md = checkmessage(md,['file provided in ''' fieldname ''': ''' field ''' does not exist']);
	end
end

%Check forcings (size and times)
if getfieldvalue(options,'timeseries',0),
	if (size(field,1)==md.mesh.numberofvertices | size(field,1)==md.mesh.numberofelements),
		if ~size(field,2)==1,
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' should have only one column as there are md.mesh.numberofvertices lines']));
		end
	elseif (size(field,1)==md.mesh.numberofvertices+1 | size(field,1)==md.mesh.numberofelements+1),
		if any(field(end,:)~=sort(field(end,:))),
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' columns should be sorted chronologically']));
		end
		if any(field(end,1:end-1)==field(end,2:end)),
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' columns must not contain duplicate timesteps']));
		end
	else
		md = checkmessage(md,getfieldvalue(options,'message',...
			['field ''' fieldname ''' should have md.mesh.numberofvertices or md.mesh.numberofvertices+1 lines']));
	end
end

%Check single value forcings (size and times)
if getfieldvalue(options,'singletimeseries',0),
	if size(field,1)==2
		if any(field(end,:)~=sort(field(end,:))),
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' columns should be sorted chronologically']));
		end
		if any(field(end,1:end-1)==field(end,2:end)),
			md = checkmessage(md,getfieldvalue(options,'message',...
				['field ''' fieldname ''' columns must not contain duplicate timesteps']));
		end
	else
		md = checkmessage(md,getfieldvalue(options,'message',...
			['field ''' fieldname ''' should have 2 lines']));
	end
end
