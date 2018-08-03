function segments=findsegments(md,varargin)
%FINDSEGMENTS - build segments model field
%
%   Optional inputs:
%      'mesh.elementconnectivity'
%
%   Usage:
%      segments=findsegments(md,varargin);

%get options
options=pairoptions(varargin{:});

%Get connectivity
mesh.elementconnectivity=getfieldvalue(options,'mesh.elementconnectivity',md.mesh.elementconnectivity);

%Now, build the connectivity tables for this mesh if not correclty done
if size(md.mesh.elementconnectivity,1)~=md.mesh.numberofelements,
	if exist(options,'mesh.elementconnectivity'),
		error(' ''mesh.elementconnectivity'' option does not have thge right size.');
	else
		mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
	end
end

%Recreate the segments
elementonboundary=double(mesh.elementconnectivity(:,3)==0);
pos=find(elementonboundary);
num_segments=length(pos);
segments=zeros(num_segments,3);
count=1;

%loop over the segments
for i=1:num_segments,

	%get current element on boundary
	el1=pos(i);

	%get elements connected to el1
	els2=mesh.elementconnectivity(el1,find(mesh.elementconnectivity(el1,:)));

	%el1 is connected to 2 other elements
	if length(els2)>1,

		%get nodes of el1
		nods1=md.mesh.elements(el1,:);

		%find the common vertices to the two elements connected to el1 (1 or 2)
		flag=intersect(md.mesh.elements(els2(1),:),md.mesh.elements(els2(2),:));

		%get the vertices on the boundary and build segment
		nods1(find(ismember(nods1,flag)))=[];
		segments(count,:)=[nods1 el1];

		%swap segment nodes if necessary
		ord1=find(nods1(1)==md.mesh.elements(el1,:));
		ord2=find(nods1(2)==md.mesh.elements(el1,:));
		if ( (ord1==1 & ord2==2) | (ord1==2 & ord2==3) | (ord1==3 & ord2==1) ),
			temp=segments(count,1);
			segments(count,1)=segments(count,2);
			segments(count,2)=temp;
		end
		segments(count,1:2)=fliplr(segments(count,1:2));
		count=count+1;

	%el1 is connected to only one element
	else
		%get nodes of el1
		nods1=md.mesh.elements(el1,:);

		%find the vertex  the el1 to not share with els2
		flag=setdiff(nods1,md.mesh.elements(els2,:));

		for j=1:3,
			nods=nods1; nods(j)=[];
			if any(ismember(flag,nods)),

				segments(count,:)=[nods el1];

				%swap segment nodes if necessary
				ord1=find(nods(1)==md.mesh.elements(el1,:));
				ord2=find(nods(2)==md.mesh.elements(el1,:));
				if ( (ord1==1 & ord2==2) | (ord1==2 & ord2==3) | (ord1==3 & ord2==1) ),
					temp=segments(count,1);
					segments(count,1)=segments(count,2);
					segments(count,2)=temp;
				end
				segments(count,1:2)=fliplr(segments(count,1:2));
				count=count+1;
			end
		end
	end
end
