function enveloppeVTK(filename,model,varargin)
% vtk export
% function enveloppeVTK(filename,model)
% creates a directory with the vtk files for displays in paraview
% only export the enveloppe result (surface and base) on trias
%
% input: filename   destination 
%                   (string)
%------------------------------------------------------------------
%        model      this is md 
%------------------------------------------------------------------
% By default only the results are exported, you can add whichever
% field you need as a string:
% add 'geometry' to export md.geometry
%
% Basile de Fleurian:

[path,name,ext]=fileparts(filename);
separator=filesep;
mkdir(filename);
IsEnveloppe=find(model.mesh.vertexonbase | model.mesh.vertexonsurface);

%get the element related variables
if dimension(model.mesh)==2,
	points=[model.mesh.x model.mesh.y zeros(model.mesh.numberofvertices,1)];
	[num_of_elt]=size(model.mesh.elements,1);
else
	points=[model.mesh.x(IsEnveloppe) model.mesh.y(IsEnveloppe) model.mesh.z(IsEnveloppe)];
	[num_of_elt]=size(find(isnan(model.mesh.lowerelements)),1)+size(find(isnan(model.mesh.upperelements)),1);
	[low_elt_num]=size(find(isnan(model.mesh.lowerelements)),1);
	[top_elt_num]=size(find(isnan(model.mesh.upperelements)),1);
end

celltype=5; %triangles
[num_of_points,dim]=size(points);
[point_per_elt]=size(model.mesh.elements,2);
tot_points=model.mesh.numberofvertices;

%this is the result structure
res_struct=model.results;
%checking for results
if (length(fields(res_struct))>0);
	%Getting all the solutions of the model
	solnames=fields(res_struct);
	num_of_sols=length(solnames);
	num_of_timesteps=1;
	%building solution structure 
	for i=1:num_of_sols
		sol_struct{i}=res_struct.(solnames{i});
		%looking for multiple time steps
		if(size(sol_struct{i},2)>num_of_timesteps);
			num_of_timesteps=size(sol_struct{i},2);
			outstep=model.timestepping.time_step*model.settings.output_frequency
    end
  end
else
	num_of_timesteps=1;
end
for step=1:num_of_timesteps;
	
	timestep=step;
	fid = fopen(strcat(path,filesep,name,filesep,'timestep.vtk',int2str(timestep),'.vtk'),'w+');
	fprintf(fid,'# vtk DataFile Version 2.0 \n');
	fprintf(fid,'Data for run %s \n',model.miscellaneous.name);
	fprintf(fid,'ASCII \n');
	fprintf(fid,'DATASET UNSTRUCTURED_GRID \n');
	
	fprintf(fid,'POINTS %d float\n',num_of_points);
	if(dim==3);
		s='%f %f %f \n';
	elseif(dim==2);
		s='%f %f \n';
  end
	P=[points zeros(num_of_points,3-dim)];
	fprintf(fid,s,P');
	
	fprintf(fid,'CELLS %d %d\n',num_of_elt,num_of_elt*(3+1));
	s='%d';
	for j=1:3
		s=horzcat(s,{' %d'});
  end
	s=cell2mat(horzcat(s,{'\n'}));

	%build the connection matrix for the top and bottom elements
	if exist('low_elt_num')
		triaconnect=zeros(num_of_elt,3);
		triaconnect(1:low_elt_num,:)=model.mesh.elements(find(isnan(model.mesh.lowerelements)),1:3);
		upshift=-min(min(model.mesh.elements(find(isnan(model.mesh.upperelements)),4:6)))+1+max(max(model.mesh.elements(find(isnan(model.mesh.lowerelements)),1:3)));
		triaconnect(1+low_elt_num:num_of_elt,:)=model.mesh.elements(find(isnan(model.mesh.upperelements)),4:6)+upshift;
		fprintf(fid,s,[(3)*ones(num_of_elt,1) triaconnect-1]');
	else
		fprintf(fid,s,[(point_per_elt)*ones(num_of_elt,1)	model.mesh.elements-1]');
  end

	fprintf(fid,'CELL_TYPES %d\n',num_of_elt);
	s='%d\n';
	fprintf(fid,s,celltype*ones(num_of_elt,1));
	fprintf(fid,'POINT_DATA %s \n',num2str(num_of_points));

	%loop over the different solution structures
	if (exist('num_of_sols'));
		for j=1:num_of_sols
			%dealing with results on different timesteps
			if(size(sol_struct{j},2)>timestep);
				timestep = step;
			else
				timestep = size(sol_struct{j},2);
	    end
			%getting the number of fields in the solution
			resfields=fields(sol_struct{j}(timestep));
			num_of_fields=length(resfields);
			%check which field is a real result and print
			for k=1:num_of_fields
				if ((numel(sol_struct{j}(timestep).(resfields{k})))==tot_points);
					%paraview does not like NaN, replacing
					nanval=find(isnan(sol_struct{j}(timestep).(resfields{k})));
					sol_struct{j}(timestep).(resfields{k})(nanval)=-9999;
					%also checking for verry small value that mess up
					smallval=(abs(sol_struct{j}(timestep).(resfields{k}))<1.0e-20);
					sol_struct{j}(timestep).(resfields{k})(smallval)=0.0;
					fprintf(fid,'SCALARS %s float 1 \n',resfields{k});
					fprintf(fid,'LOOKUP_TABLE default\n');
					s='%e\n';
					fprintf(fid,s,sol_struct{j}(timestep).(resfields{k})(IsEnveloppe));
		    end		
	    end 
	  end
  end
	%loop on arguments, if something other than result is asked, do
	%it now
	for j= 1:nargin-2
		res_struct=model.(varargin{j});
		fieldnames=fields(res_struct);
		num_of_fields=length(fieldnames);
		for k=1:num_of_fields
			if ((numel(res_struct.(fieldnames{k})))==tot_points);
				%paraview does not like NaN, replacing
				nanval=find(isnan(res_struct.(fieldnames{k})));
				res_struct.(fieldnames{k})(nanval)=-9999;
				%also checking for verry small value that mess up
				smallval=(abs(res_struct.(fieldnames{k}))<1.0e-20);
				res_struct.(fieldnames{k})(smallval)=0.0;
				fprintf(fid,'SCALARS %s float 1 \n',fieldnames{k});
				fprintf(fid,'LOOKUP_TABLE default\n');
				s='%e\n';
				fprintf(fid,s,res_struct.(fieldnames{k})(IsEnveloppe));
				%check for forcings	
			elseif (size(res_struct.(fieldnames{k}),1)==tot_points+1);
				%paraview does not like NaN, replacing
				nanval=find(isnan(res_struct.(fieldnames{k})));
				res_struct.(fieldnames{k})(nanval)=-9999;
				%also checking for verry small value that mess up
				smallval=(abs(res_struct.(fieldnames{k}))<1.0e-20);
				res_struct.(fieldnames{k})(smallval)=0.0;
				if (size(res_struct.(fieldnames{k}),2)==num_of_timesteps),
					fprintf(fid,'SCALARS %s float 1 \n',fieldnames{k});
					fprintf(fid,'LOOKUP_TABLE default\n');
					s='%e\n';
					fprintf(fid,s,res_struct.(fieldnames{k})(IsEnveloppe,timestep));
				else,
					%forcing and results not on the same timestep,need some treatment
					fprintf(fid,'SCALARS %s float 1 \n',fieldnames{k});
					fprintf(fid,'LOOKUP_TABLE default\n');
					index=1
					currenttime=((timestep-1)*outstep)+model.timestepping.start_time+model.timestepping.time_step
					while (res_struct.(fieldnames{k})(end,index)<=currenttime);
						index=index+1
		      end
					uptime=res_struct.(fieldnames{k})(end,index);
					uplim=res_struct.(fieldnames{k})(IsEnveloppe,index);
					uptime
					while (res_struct.(fieldnames{k})(end,index)>=currenttime);
						index=index-1
		      end
					lowtime=res_struct.(fieldnames{k})(end,index);
					lowlim=res_struct.(fieldnames{k})(IsEnveloppe,index);
					lowtime
					interp=lowlim+(uplim-lowlim)*((currenttime-lowtime)/(uptime-lowtime))
					s='%e\n';
					fprintf(fid,s,interp);
				end	
		  end		
		end 
	end
	fclose(fid);
end
