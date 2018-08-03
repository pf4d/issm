function dvar=QmuSetupVariables(md,dvar,variables)

%get descriptor
descriptor=variables.descriptor;

%decide whether this is a distributed variable, which will drive whether we expand it into npart values,
%or if we just carry it forward as is. 

%ok, key off according to type of descriptor:
if strncmp(descriptor,'scaled_',7),
	%we have a scaled variable, expand it over the partition.

	if isa(variables,'uniform_uncertain'),
		if (length(variables.lower)>md.qmu.numberofpartitions || length(variables.upper)>md.qmu.numberofpartitions)
			error('QmuSetupDesign error message: stddev should be either a scalar or a ''npart'' length vector');
		end
	elseif isa(variables,'normal_uncertain'),
		if length(variables.stddev)>md.qmu.numberofpartitions,
			error('QmuSetupDesign error message: stddev should be either a scalar or a ''npart'' length vector');
		end
	end

	%ok, dealing with semi-discrete distributed variable. Distribute according to how many 
	%partitions we want

	for j=1:md.qmu.numberofpartitions
		dvar(end+1)           =variables;
		dvar(end  ).descriptor=sprintf('%s_%d',variables.descriptor,j);
		if isa(variables,'uniform_uncertain'),
			if length(variables.lower)>1,
				dvar(end  ).lower=variables.lower(j);
			end
			if length(variables.upper)>1,
				dvar(end  ).upper=variables.upper(j);
			end
		elseif isa(variables,'normal_uncertain'),
			if length(variables.stddev)>1,
				dvar(end  ).stddev=variables.stddev(j);
			end
		end
	end

else
	dvar(end+1)=variables;
end
