function plot_transient_movie(md,options,width,i)
%PLOT_TRANSIENT_MOVIE - plot a transient result as a movie
%   Usage:
%      plot_transient_movie(md,options,width,i);
%
%   See also: PLOTMODEL, PLOT_UNIT, PLOT_MANAGER

	%prepare subplot
	subplot(width,width,i); 

	%xlim
	if exist(options,'transient_movie_field'),
		field=getfieldvalue(options,'transient_movie_field');
	else
		disp('List of available fields:');
		F=fields(md.results.TransientSolution(1));
		num = [];
		for i=1:numel(F),
			if ~strcmp(F{i},'time') & ...
				~strcmp(F{i},'step') & ...
				~strcmp(F{i},'errlog') & ...
				~strcmp(F{i},'outlog') & ...
				~strcmp(F{i},'MaxIterationConvergenceFlag') & ...
				~strcmp(F{i},'SolutionType'),
				disp(['   ' num2str(i) ': ' F{i} ]);
				num = [num i];
			end
		end
		choice=input(['please enter the field number? (between ' num2str(min(num)) ' and ' num2str(max(num)) ')  ']);
		field =  F{choice};
	end

	results=md.results.TransientSolution;
	%loop over the time steps
	if exist(options,'transient_movie_limit'),
		limit=getfieldvalue(options,'transient_movie_limit');
		steps=[limit(1):limit(end)];
	else
		steps=1:length(results);
	end

	%calculate caxis
	if ~exist(options,'caxis'),
		range = [Inf -Inf];
		for i=steps
			[data datatype]=processdata(md,results(i).(field),options);
			range(1) = min(range(1),min(data));
			range(2) = max(range(2),max(data));
		end
		options=addfielddefault(options,'caxis',range);
	end

	%display movie
	nstep=1;
	deltat = getfieldvalue(options,'pause',.5);
	for i=steps

		if ~isempty(results(i).(field)),
			%process data
			[x y z elements is2d isplanet]=processmesh(md,results(i).(field),options);
			[data datatype]=processdata(md,results(i).(field),options);

			clf;
			titlestring=[field ' at time ' num2str(results(i).time) ' year'];
			plot_unit(x,y,z,elements,data,is2d,isplanet,datatype,options)
			apply_options_movie(md,options,titlestring);

			if exist(options,'transient_movie_output'),
				set(gcf,'Renderer','zbuffer','color','white'); %fixes a bug on Mac OS X (not needed in future Matlab version)
				if nstep==1,
					%initialize images and frame
					frame=getframe(gcf);
					[images,map]=rgb2ind(frame.cdata,256,'nodither');
					images(1,1,1,length(steps))=0;
				else
					frame=getframe(gcf);
					images(:,:,1,nstep) = rgb2ind(frame.cdata,map,'nodither');
				end
			else
				pause(deltat)
			end
			nstep=nstep+1;
		end
	end

	%output movie if requested.
	if exist(options,'transient_movie_output'),
		filename=getfieldvalue(options,'transient_movie_output');
		imwrite(images,map,filename,'DelayTime',getfieldvalue(options,'transient_movie_time',2),'LoopCount',inf)
	end

end %function

function apply_options_movie(md,options,titlestring)
	%apply options
	options=changefieldvalue(options,'title',titlestring);
	options=addfielddefault(options,'colorbar',1);
	applyoptions(md,[],options);
end
