%ga: grounded area
%iv: ice volume
%ivaf: ice volume above floatation
%GLy40 : grounding line position @ y=40km
%nelem : number of elements

function [ga iv ivaf GLy40 nelem t] = ice_evolution(md),

	ga			= [];
	iv			= [];
	ivaf		= [];
	GLy40		= [];
	nelem		= [];
	t			= [];
	nsteps	= length(md.results.TransientSolution);

	for i=1:nsteps,
		ga(i)			= md.results.TransientSolution(i).GroundedArea;
		iv(i)			= md.results.TransientSolution(i).IceVolume;
		ivaf(i)		= md.results.TransientSolution(i).IceVolumeAboveFloatation;
		nelem(i)		= size(md.results.TransientSolution(i).MeshElements,1);
		t(i)			= md.results.TransientSolution(i).time;	
		%find GL position at y=40km
		[glx gly]	= gl_position(md,i,0);
		pos=find(gly<45000 & gly > 35000);
		x=gly(pos);
		v=glx(pos);
		xq=[38000:100:42000];
		vq = interp1(x,v,xq,'linear');
		pos=find(xq==40000);
		GLy40(i)=vq(pos);
	end

end
