
fid=fopen('runme.m','r');

tline = fgets(fid);
count=1;
while ischar(tline)
	tline = fgets(fid);
	if length(tline)>16,
		if strcmpi(tline(1:16),'if perform(org,'''),
			disp(sprintf('%i: %s',count,tline(17:end-4)));
			count=count+1;
		end
	end
end
