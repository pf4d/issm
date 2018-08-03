function ids=IdFromString(string),
%IDFROMSTRING - output ids from a given string
%
%   Usage:
%      ids=IdFromString(string);
%
%   Examples:
%      ids=IdFromString('Parallel');
%      ids=IdFromString('79North');
%      ids=IdFromString('*');          %Print all tests

%Check input
if ~ischar(string)
	error('IdFromString error message: input argument is not a string');
end

%Initialize output
ids=[];

%Grep string
if strcmp(string,'*'),
	[dummy ids_raw]=system(['find ./ -name "test[0-9]*.m" | xargs grep "%Test Name:" | sed -e "s/test/ /g" -e "s/\.m:/ /g" | awk {''print $2''}']);
else
	[dummy ids_raw]=system(['find ./ -name "test[0-9]*.m" | xargs grep "%Test Name:" | grep ' string ' | sed -e "s/test/ /g" -e "s/\.m:/ /g" | awk {''print $2''}']);
end

%return if no test found
if isempty(ids_raw),
	disp(['No test matches ''' string '''' ]);
	return
end

%Process string (delete return carriage);
ids_raw=strsplit_strict(ids_raw,char(10));
ids_raw=ids_raw(1:end-1);
for i=1:length(ids_raw),
	eval(['ids=[ids ' ids_raw{i} '];']); 
end
ids=sort(ids);

%Display names
disp([ num2str(length(ids)) ' tests match ''' string '''']);
for i=1:length(ids)
	disp([ '   ' num2str(ids(i)) ' : ' IdToName(ids(i)) ]);
end
