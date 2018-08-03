steps=[1];

if any(steps==1) %Transient Run #1

	md = loadmodel('../Pig/Models/PIG_Control_drag');	

	md.inversion.iscontrol=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
	md.transient.ismovingfront=0;
	md.transient.isthermal=0;
	
	pos=find(md.mask.groundedice_levelset<0);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=25*ones(md.mesh.numberofvertices,1);
	
	md.timestepping.time_step=0.1;
	md.timestepping.final_time=10;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

	md=solve(md,'Transient');

	% Save model
	save ./Models/PIG_Transient md;
end

if any(steps==2) %High Melt #2
	md = loadmodel('./Models/PIG_Transient');	

	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=60*ones(md.mesh.numberofvertices,1);
	
	md.timestepping.time_step=0.1;
	md.timestepping.final_time=10;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

	md=solve(md,'Transient');

	save ./Models/PIG_HighMelt md;
end

if any(steps==3) %Ice Front retreat
	md = loadmodel('./Models/PIG_Transient');	

	md2=extract(md,'~FrontRetreat.exp');

	md2=SetMarineIceSheetBC(md2);

	md2.basalforcings.groundedice_melting_rate=zeros(md2.mesh.numberofvertices,1);
	md2.basalforcings.floatingice_melting_rate=25*ones(md2.mesh.numberofvertices,1);

	md2.timestepping.time_step=0.1;
	md2.timestepping.final_time=10;
	md2.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

	md2=solve(md2,'Transient');

	save ./Models/PIG_FrontRetreat md2;
end

if any(steps==4) %High surface mass balance #3
	%Load model

	%Change external forcing basal melting rate and surface mass balance)
	
	%Refine time steps and time span of the simulation

	%Request additional outputs

	%Solve

	%Save model

end
