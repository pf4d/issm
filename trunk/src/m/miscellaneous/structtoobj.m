function obj=structtoobj(obj,S),
%Convert struct to object

	%Get object and structure fields
	structfields=fields(S);
	objprops    =properties(class(obj));

	%recover object properties
	for i=1:length(structfields),
		fieldname =structfields{i};
		if ismember(fieldname,objprops),
			fieldvalue=getfield(S,fieldname);
			obj=setfield(obj,fieldname,fieldvalue);
		end
	end
end
