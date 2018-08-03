function export_netCDF(md,filename)	
	
%Now going on Real treatment
	if exist(filename),
		disp(sprintf('File %s allready exist', filename));
		prompt = 'Give a new name or "delete" to replace: ';
		newname = input(prompt,'s');
		if strcmp(newname,'delete')
			delete(filename)
		else
			disp(sprintf('New file name is %s ', newname));
			filename=newname
	  end
  end
	%open file and write description
	mode = netcdf.getConstant('NC_NETCDF4');
	mode = bitor(mode,netcdf.getConstant('NC_NOCLOBBER'));%NOCLOBBER to avoid overwrite
	ncid = netcdf.create(filename,mode);
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title',['Results for run ' md.miscellaneous.name]);
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date',['Created ' datestr(now)]);
	
	%gather geometry and timestepping as dimensions
	resfields=fieldnames(md.results);
	Duration=size(eval(['md.results. ' resfields{1} ]),2);
	if Duration>0,
		StepNum=Duration;
	else
		StepNum=1;
  end							

   dimlist=[40,2,md.mesh.numberofelements,md.mesh.numberofvertices,size(md.mesh.elements,2)];
 
	%define netcdf dimensions
	DimSize(1).index=netcdf.defDim(ncid,'Dimension1',StepNum);
	[DimSize(1).name,DimSize(1).value]=netcdf.inqDim(ncid,DimSize(1).index);
	DimValue(1)=DimSize(1).value;
	for i=1:5
		if sum(dimlist(i) == DimValue) == 0
			DimSize(i+1).index=netcdf.defDim(ncid,['Dimension' num2str(i+1)],dimlist(i));
			[DimSize(i+1).name,DimSize(i+1).value]=netcdf.inqDim(ncid,DimSize(i+1).index);
			DimValue(i+1)=DimSize(i+1).value;
		end
	end

	typelist=[{'numeric'} {'logical'} {'string'} {'char'} {'cell'}];
 
	%get all model classes and create respective groups
	groups=fieldnames(md);
	for i=1:length(groups),
		disp(sprintf('group name in tree %s ',groups{i}));
		groupID=netcdf.defGrp(ncid,groups{i});
		%In each group gather the fields of the class
		groupfields=fields(md.(groups{i}));
		for j=1:length(groupfields),
			Var=md.(groups{i}).(groupfields{j});
			if isa(Var,'cell')
				Stdlist=false;
				if length(Var) == 0
					Stdlist=true;
				else
					for k=1:length(typelist)
						if isa(Var{1},typelist{k})
							Stdlist=true;
						end
					end
				end

				netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),'classtype',class(md.(groups{i})));
				if(Stdlist)
					disp(sprintf('=====Field name in tree %s ',groupfields{j}));
					[DimSize,DimValue]=DefCreateVar(ncid,Var,groupID,groupfields{j},DimSize,DimValue);
				else
					listsize=length(Var);
					subgroupID=netcdf.defGrp(groupID,groupfields{j});
					netcdf.putAtt(subgroupID,netcdf.getConstant('NC_GLOBAL'),'classtype',class(Var));
					for l=1:listsize
						if isprop(Var{l},'name')
							lname=Var{l}.name;
						elseif isprop(Var{l},'step')
							lname=Var{l}.step
						else 
							lname=[class(Var{l}) int2str(l)];
						end
						listgroupID=netcdf.defGrp(subgroupID,lname);
						netcdf.putAtt(listgroupID,netcdf.getConstant('NC_GLOBAL'),'classtype',class(Var{l}));
						subfields=fields(Var{l});
						for m=1:length(subfields)
							if ~strcmp(subfields{m},'outlog')
								[DimSize,DimValue]=DefCreateVar(ncid,Var{l}.(subfields{m}),listgroupID,subfields{m},DimSize,DimValue);
							end
						end
					end
				end
			elseif isa(Var,'struct') && ~strcmp(groupfields{j},'bamg')
				classtype=class(md.(groups{i}));
				if strcmp(classtype,'struct')
					classtype=groups{i};
				end
				netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),'classtype',classtype);
				if length(Var)>1
					listsize=length(Var);
					subgroupID=netcdf.defGrp(groupID,groupfields{j});
					classtype=class(Var);
					if strcmp(classtype,'struct')
						classtype=groups{i};
					end
					netcdf.putAtt(subgroupID,netcdf.getConstant('NC_GLOBAL'),'classtype',classtype);
					for l=1:listsize
						if isfield(Var(l),'step')
							lname=[int2str(Var(l).step)];
						else
							lname=[class(Var(l)) int2str(l)];
						end
						classtype=class(Var(l));
						if strcmp(classtype,'struct')
							classtype=groups{i};
						end
						listgroupID=netcdf.defGrp(subgroupID,lname);
						netcdf.putAtt(listgroupID,netcdf.getConstant('NC_GLOBAL'),'classtype',classtype);
						subfields=fields(Var(l));
						for m=1:length(subfields)
							if ~strcmp(subfields{m},'outlog')
								[DimSize,DimValue]=DefCreateVar(ncid,Var(l).(subfields{m}),listgroupID,subfields{m},DimSize,DimValue);
							end
						end
					end
				else
					subgroupID=netcdf.defGrp(groupID,groupfields{j});
					classtype=class(Var);
					if strcmp(classtype,'struct')
						classtype=groups{i};
					end
					netcdf.putAtt(subgroupID,netcdf.getConstant('NC_GLOBAL'),'classtype',classtype);
					subfields=fields(Var);
					for m=1:length(subfields)
						if ~strcmp(subfields{m},'outlog')
							[DimSize,DimValue]=DefCreateVar(ncid,Var.(subfields{m}),subgroupID,subfields{m},DimSize,DimValue);
						end
					end
				end
			else
				netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),'classtype',class(md.(groups{i})));
				[DimSize,DimValue]=DefCreateVar(ncid,Var,groupID,groupfields{j},DimSize,DimValue);
			end
		end
 end
 netcdf.close(ncid);
end

function [DimSize,DimValue]=DefCreateVar(ncid,Var,groupID,field,DimSize,DimValue,last,md,midfield)
	varclass=class(Var);
	varsize=size(Var);
	varlength=length(Var);
	if isa(Var,'logical'),
		if Var,
			LogicString='True';
		else,
			LogicString='False';
  	end
		netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),field,LogicString);
	elseif isa(Var,'char'),
		netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),field,Var);
	elseif isa(Var,'double'), %dealing with arrays
		[dims,DimSize,DimValue]=GetDims(ncid,Var,DimSize,DimValue);
 		varid = netcdf.defVar(groupID,field,'NC_DOUBLE',dims);
		if length(Var)==0,
			netcdf.putVar(groupID,varid,NaN);
		else
			netcdf.putVar(groupID,varid,Var);
		end
	elseif isa(Var,'cell'),
		[dims,DimSize,DimValue]=GetDims(ncid,Var,DimSize,DimValue);
		%dirty hack to be able to pass strings
		varid = netcdf.defVar(groupID,field,'NC_CHAR',dims);
		if length(Var)==0,
			netcdf.putVar(groupID,varid,0,9,'emptycell')
		else
			for i=1:length(Var),
				if length(Var)>1,
					endpoint=[1,min(length(Var{i}),40)];
					startpoint=[1 0];
				else
					endpoint=min(length(Var{i}),40);
					startpoint=0;
				end
				if length(Var{i})>40,
					netcdf.putVar(groupID,varid,startpoint,extent,Var{i}(1:40))
					disp(sprintf('some variable have been truncated'));
				else
					netcdf.putVar(groupID,varid,startpoint,endpoint,Var{i})
				end
			end
		end
	elseif isa(Var,'struct'),
		%Start by getting the structure fields and size
		locfields=fields(Var);
		[dims,DimSize,DimValue]=GetDims(ncid,Var,DimSize,DimValue);
		varid = netcdf.defVar(groupID,field,'NC_CHAR',dims);
		if length(locfields)==0,
			netcdf.putVar(groupID,varid,[0,0],[1,11],'emptystruct')
		else
			for i=1:length(locfields),
				for j=1:2,
					if j==1,
						CharVar=locfields{i};
						if length(CharVar)==0
							CharVar='emptystruct';
						end
						startpoint=[i-1,0,0];
					else
						if isa(Var.(locfields{i}),'char'),
							CharVar=Var.(locfields{i});
						else
							CharVar=num2str(Var.(locfields{i}));
						end
						if length(CharVar)==0
							CharVar='emptystruct';
						end
						startpoint=[i-1,1,0];
					end

					extent=[1,1,min(length(CharVar),40)];
					if length(CharVar)>40,
						netcdf.putVar(groupID,varid,startpoint,extent,CharVar(1:40))
						disp(sprintf('some variable have been truncated'));
					else
						netcdf.putVar(groupID,varid,startpoint,extent,CharVar)
					end
				end
			end
		end
	else
		disp(sprintf('no support for class %s of field %s',varclass,field));
  end
	return
end

function [dims,DimSize,DimValue]=GetDims(ncid,Var,DimSize,DimValue)
	dims=[];
	if isa(Var,'cell'),
		varsize=size(Var');
	elseif isa(Var,'struct')
		varsize=length(fields(Var));
	else
		varsize=size(Var);
	end
	dim=sum(varsize>1);
	if dim>0
		for i=1:dim
			indsize=find(varsize(i)==DimValue);
			if length(indsize)>0
				dims=[dims DimSize(indsize).index];
			else
				indsize=length(DimSize)+1;
				DimSize(indsize).index=netcdf.defDim(ncid,['Dimension' num2str(indsize)],varsize(i));
				[DimSize(indsize).name,DimSize(indsize).value]=netcdf.inqDim(ncid,DimSize(indsize).index);
				DimValue(indsize)=DimSize(indsize).value;
				dims=[dims DimSize(indsize).index];
			end
		end
	end
	%if we have a cell variable we need to add a stringlength dimension 
	if isa(Var,'struct'),
		if DimValue(3)~=2
			if DimValue(2)~=2
				dims=[dims DimSize(1).index];
			else
				dims=[dims DimSize(2).index];
			end
		else
			dims=[dims DimSize(3).index];
		end
	end
	if isa(Var,'cell') || isa(Var,'struct'),
		if DimValue(2)~=40
			dims=[dims DimSize(1).index];
		else
			dims=[dims DimSize(2).index];
		end
	end
end
