function varargout=runme(varargin)
%RUNME - test deck for ISSM nightly runs
%
%   In a test deck directory (tests/Vertification/NightlyRun for example)
%   The following command will launch all the existing tests:
%   >> runme
%   To run the tests 101 and 102:
%   >> runme('id',[101 102])
%   etc...
%
%   Available options:
%      'id'            followed by the list of ids requested
%      'exclude'       ids to be excluded from the test
%      'benchmark'     'all' (all of them)
%                      'nightly' (nightly run)
%                      'ismip'  : validation of ismip-hom tests
%                      'eismint': validation of eismint tests
%                      'thermal': validation of thermal tests
%                      'mesh'   : validation of mesh tests
%                      'adolc'   : validation of adolc tests
%                      'slr'   : validation of slr tests
%                      'qmu'   : validation of dakota tests
%                      ...
%      'procedure'     'check' :   run the test (default)
%                      'update':   update the archive
%                      'valgrind': check for memory leaks (default value of md.debug.valgrind needs to be changed manually)
%      'stoponerror'   1 or 0
%
%   Usage:
%      runme(varargin);
%
%   Examples:
%      runme;
%      runme('exclude',101);
%      runme('id',102,'procedure','update');

%Check inputs
% {{{
if nargout>1
	help runme
	error('runme error message: bad usage');
end

%recover options
options=pairoptions(varargin{:});
% }}}

%Process options
%GET benchmark {{{
benchmark=getfieldvalue(options,'benchmark','nightly');
if ~ismember(benchmark,{'all','nightly','ismip','eismint','thermal','mesh','validation','tranforcing','adolc','slr','qmu'})
	disp('runme warning: benchmark not supported, defaulting to test ''nightly''')
	benchmark='nightly';
end
% }}}
%GET procedure {{{
procedure=getfieldvalue(options,'procedure','check');
if ~ismember(procedure,{'check','update','valgrind'})
	disp('runme warning: procedure not supported, defaulting to test ''check''')
	procedure='check';
end
% }}}
%GET output {{{
output=getfieldvalue(options,'output','none');
if ~ismember(output,{'nightly','none'})
	disp('runme warning: output not supported, defaulting to test ''none''')
	output='none';
end
% }}}
%GET RANK and NUMPROCS for multithreaded runs  {{{
rank=getfieldvalue(options,'rank',1);
numprocs=getfieldvalue(options,'numprocs',1);
if (numprocs<rank), numprocs=1; end
% }}}
%GET ids  {{{
flist=dir;%use dir, as it seems to act OS independent
list_ids=[];
for i=1:numel(flist),
	if ( strncmp(flist(i).name,'test',4) &...                         %File name must start with 'test'
			strncmp(fliplr(flist(i).name),fliplr('.m'),2)&...           %File name must end by '.m'
			~strcmp(flist(i).name,'test.m'))                            %File name must be different than 'test.m'
		id=str2num(flist(i).name(5:end-2));
		if isempty(id),
			disp(['WARNING: ignore file ' flist(i).name ]);
		else
			list_ids(end+1)=eval(flist(i).name(5:end-2));                  %Keep test id only (skip 'test' and '.m')
		end
	end
end
[i1,i2]=parallelrange(rank,numprocs,length(list_ids));               %Get tests for this cpu only
list_ids=list_ids(i1:i2);

test_ids=getfieldvalue(options,'id',list_ids);
test_ids=intersect(test_ids,list_ids);
% }}}
%GET exclude {{{
exclude_ids=getfieldvalue(options,'exclude',[]);
exclude_ids=[exclude_ids];
pos=find(ismember(test_ids,exclude_ids));
test_ids(pos)=[];
% }}}
%Process Ids according to benchmarks{{{
if strcmpi(benchmark,'nightly'),
	test_ids=intersect(test_ids,[1:999]);
elseif strcmpi(benchmark,'validation'),
	test_ids=intersect(test_ids,[1001:1999]);
elseif strcmpi(benchmark,'ismip'),
	test_ids=intersect(test_ids,[1101:1199]);
elseif strcmpi(benchmark,'eismint'),
	test_ids=intersect(test_ids,[1201:1299]);
elseif strcmpi(benchmark,'thermal'),
	test_ids=intersect(test_ids,[1301:1399]);
elseif strcmpi(benchmark,'mesh'),
	test_ids=intersect(test_ids,[1401:1499]);
elseif strcmpi(benchmark,'tranforcing'),
	test_ids=intersect(test_ids,[1501:1502]);
elseif strcmpi(benchmark,'referential'),
	test_ids=intersect(test_ids,[1601:1602]);
elseif strcmpi(benchmark,'slr'),
	test_ids=intersect(test_ids,[2001:2500]);
elseif strcmpi(benchmark,'adolc'),
	test_ids=intersect(test_ids,[3001:3200]);
elseif strcmpi(benchmark,'qmu'),
	test_ids=intersect(test_ids,[218 234 235 412:414 417 418 420]);
end
% }}}

%Loop over tests and launch sequence
root=pwd;
for id=test_ids,
	disp(sprintf('%s%i%s','----------------starting:',id,'-----------------------'));
	try,
		%Execute test
		eval(['cd ' root ]);
		id_string='N/A';
		id_string=IdToName(id);
		eval(['test' num2str(id)]);

		%UPDATE ARCHIVE?
		archive_name=['Archive' num2str(id) ];
		if strcmpi(procedure,'update'),
			delete(['../Archives/' archive_name '.arch'])
			for k=1:length(field_names),
				field=field_values{k};
				archwrite(['../Archives/' archive_name '.arch'],[archive_name '_field' num2str(k)], field);
			end
			disp(sprintf(['File ./../Archives/' archive_name '.arch saved\n']));

		%CHECK for memory leaks?
		elseif strcmpi(procedure,'valgrind'),
			fields = fieldnames(md.results);
			for i=1:numel(fields)
				if ~isfield(md.results.(fields{i}),'errlog'),
					disp(['Skipping ' fields{i}]);
					continue;
				else
					disp(['Extracting results of ' fields{i}]);
				end
				results = md.results.(fields{i});
				errlog  = cellstr(results(1).errlog);

				%Check leaks
				lines  = strfind(errlog,'definitely lost:');
				lines  = find(~cellfun(@isempty,lines));
				leaks   = 0;
				for j=1:numel(lines)
					Line    = errlog(lines(j));
					Numbers = sscanf(Line{1},'==%i==   definitely lost: %s bytes in %i blocks',[1 Inf]);
					leaks   = leaks + str2num(strrep(char(Numbers(2:end-1)),',',''));
				end
				%Check conditional jumps
				lines  = strfind(errlog,'Conditional jump or move depends on uninitialised value');
				lines  = find(~cellfun(@isempty,lines));
				jumps   = numel(lines);
				%Check invalid read/write
				lines  = strfind(errlog,'Invalid');
				lines  = find(~cellfun(@isempty,lines));
				inval  = numel(lines);
				if leaks==0,
					disp(sprintf(['SUCCESS difference: 0 < 0 test id: %i test name: %s field: valgrind mem. leaks'],id,id_string));
				else
					disp(sprintf(['ERROR   difference: %i > 0 test id: %i test name: %s field: valgrind mem. leaks'],leaks,id,id_string));
					disp('STOP');
					return;
				end
				if jumps==0,
					disp(sprintf(['SUCCESS difference: 0 < 0 test id: %i test name: %s field: valgrind cond. jumps'],id,id_string));
				else
					disp(sprintf(['ERROR   difference: %i > 0 test id: %i test name: %s field: valgrind cond. jumps'],jumps,id,id_string));
					disp('STOP');
					return;
				end
				if inval==0,
					disp(sprintf(['SUCCESS difference: 0 < 0 test id: %i test name: %s field: valgrind invalid read/write'],id,id_string));
				else
					disp(sprintf(['ERROR   difference: %i > 0 test id: %i test name: %s field: valgrind invalid read/write'],inval,id,id_string));
					disp('STOP');
					return;
				end
			end

		%ELSE: CHECK TEST
		else,
			for k=1:length(field_names),

				try,
					%Get field and tolerance
					field=field_values{k};
					fieldname=field_names{k};
					tolerance=field_tolerances{k};

					%compare to archive
					%our output is in the correct order (n,1) or (1,1), so we do not need to transpose again
					archive_cell=archread(['../Archives/' archive_name '.arch'],[archive_name '_field' num2str(k)]);
					archive=archive_cell{1};
					error_diff=full(max(abs(archive(:)-field(:)))/(max(abs(archive(:)))+eps));

					%disp test result
					if (error_diff>tolerance | isnan(error_diff));
						disp(sprintf(['ERROR   difference: %-7.2g > %7.2g test id: %i test name: %s field: %s'],...
							error_diff,tolerance,id,id_string,fieldname));
						if(getfieldvalue(options,'stoponerror',0)), disp('STOP'); return; end
					else
						disp(sprintf(['SUCCESS difference: %-7.2g < %7.2g test id: %i test name: %s field: %s'],...
							error_diff,tolerance,id,id_string,fieldname));
					end

				catch me2

					%something went wrong, print failure message:
					message=getReport(me2);
					fprintf('%s',message);
					if strcmpi(output,'nightly')
						fid=fopen([issmdir() '/nightlylog/matlaberror.log'], 'at');
						fprintf(fid,'%s',message);
						fprintf(fid,'\n------------------------------------------------------------------\n');
						fclose(fid);
						disp(sprintf(['FAILURE difference: N/A test id: %i test name: %s field: %s'],id,id_string,fieldname));
					else
						disp(sprintf(['FAILURE difference: N/A test id: %i test name: %s field: %s'],id,id_string,fieldname));
						fprintf('%s',message);
						if(getfieldvalue(options,'stoponerror',0)), disp('STOP'); return; end
					end
					continue;
				end
			end
		end
	catch me,

		%something went wrong, print failure message:
		message=getReport(me);
		fprintf('%s',message);
		if strcmpi(output,'nightly')
			fid=fopen([issmdir() '/nightlylog/matlaberror.log'], 'at');
			fprintf(fid,'%s',message);
			fprintf(fid,'\n------------------------------------------------------------------\n');
			fclose(fid);
			disp(sprintf(['FAILURE difference: N/A test id: %i test name: %s field: %s'],id,id_string,'N/A'));
		else
			disp(sprintf(['FAILURE difference: N/A test id: %i test name: %s field: %s'],id,id_string,'N/A'));
			rethrow(me);
			if(getfieldvalue(options,'stoponerror',0)), disp('STOP'); return; end
		end
	end
	disp(sprintf('%s%i%s','----------------finished:',id,'-----------------------'));
end
