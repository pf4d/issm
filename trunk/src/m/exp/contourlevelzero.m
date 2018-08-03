function contours=contourlevelzero(md,mask,level)
%CONTOURLEVELZERO - figure out the zero level (or offset thereof, specified by the level value) 
%                   of a vectorial mask, and vectorialize it into an exp or shp compatible 
%structure.
%
%   Usage:
%      contours=contourlevelzero(md,mask,level)
%
%   See also: PLOT_CONTOUR

%process data 
if dimension(md.mesh)==3,
	error('contourlevelzero error message: routine not supported for 3d meshes, project on a layer');
end

x=md.mesh.x;
y=md.mesh.y;
index=md.mesh.elements;

if isempty(mask), error('mask provided is empty'); end
	if length(mask)~=md.mesh.numberofvertices, error('mask provided should be specified at the vertices of the mesh'); end

%initialization of some variables
numberofelements=size(index,1);
elementslist=1:numberofelements;
c=[];
h=[];

%get unique edges in mesh
%1: list of edges
edges=[index(:,[1,2]); index(:,[2,3]); index(:,[3,1])];
%2: find unique edges
[edges,I,J]=unique(sort(edges,2),'rows');
%3: unique edge numbers
vec=J;
%4: unique edges numbers in each triangle (2 triangles sharing the same edge will have
%   the same edge number)
edges_tria=[vec(elementslist), vec(elementslist+numberofelements), vec(elementslist+2*numberofelements)];

%segments [nodes1 nodes2]
Seg1=index(:,[1 2]);
Seg2=index(:,[2 3]);
Seg3=index(:,[3 1]);

%segment numbers [1;4;6;...]
Seg1_num=edges_tria(:,1);
Seg2_num=edges_tria(:,2);
Seg3_num=edges_tria(:,3);

%value of data on each tips of the segments
Data1=mask(Seg1);
Data2=mask(Seg2);
Data3=mask(Seg3);

%get the ranges for each segment
Range1=sort(Data1,2);
Range2=sort(Data2,2);
Range3=sort(Data3,2);

%find the segments that contain this value
pos1=(Range1(:,1)<level & Range1(:,2)>level);
pos2=(Range2(:,1)<level & Range2(:,2)>level);
pos3=(Range3(:,1)<level & Range3(:,2)>level);

%get elements
poselem12=(pos1 & pos2);
poselem13=(pos1 & pos3);
poselem23=(pos2 & pos3);
poselem=find(poselem12 | poselem13 | poselem23);
numelems=length(poselem);

%if no element has been flagged, skip to the next level
if numelems==0,
	warning('contourlevelzero warning message: no elements found with corresponding level value in mask');
	contours=struct([]);
	return;
end

%go through the elements and build the coordinates for each segment (1 by element)
x1=zeros(numelems,1);
x2=zeros(numelems,1);
y1=zeros(numelems,1);
y2=zeros(numelems,1);
edge_l=zeros(numelems,2);

for j=1:numelems,

	weight1=(level-Data1(poselem(j),1))/(Data1(poselem(j),2)-Data1(poselem(j),1));
	weight2=(level-Data2(poselem(j),1))/(Data2(poselem(j),2)-Data2(poselem(j),1));
	weight3=(level-Data3(poselem(j),1))/(Data3(poselem(j),2)-Data3(poselem(j),1));

	if poselem12(poselem(j));

		x1(j)=x(Seg1(poselem(j),1))+weight1*(x(Seg1(poselem(j),2))-x(Seg1(poselem(j),1)));
		x2(j)=x(Seg2(poselem(j),1))+weight2*(x(Seg2(poselem(j),2))-x(Seg2(poselem(j),1)));
		y1(j)=y(Seg1(poselem(j),1))+weight1*(y(Seg1(poselem(j),2))-y(Seg1(poselem(j),1)));
		y2(j)=y(Seg2(poselem(j),1))+weight2*(y(Seg2(poselem(j),2))-y(Seg2(poselem(j),1)));
		edge_l(j,1)=Seg1_num(poselem(j));
		edge_l(j,2)=Seg2_num(poselem(j));

	elseif poselem13(poselem(j)),

		x1(j)=x(Seg1(poselem(j),1))+weight1*(x(Seg1(poselem(j),2))-x(Seg1(poselem(j),1)));
		x2(j)=x(Seg3(poselem(j),1))+weight3*(x(Seg3(poselem(j),2))-x(Seg3(poselem(j),1)));
		y1(j)=y(Seg1(poselem(j),1))+weight1*(y(Seg1(poselem(j),2))-y(Seg1(poselem(j),1)));
		y2(j)=y(Seg3(poselem(j),1))+weight3*(y(Seg3(poselem(j),2))-y(Seg3(poselem(j),1)));
		edge_l(j,1)=Seg1_num(poselem(j));
		edge_l(j,2)=Seg3_num(poselem(j));

	elseif poselem23(poselem(j)),

		x1(j)=x(Seg2(poselem(j),1))+weight2*(x(Seg2(poselem(j),2))-x(Seg2(poselem(j),1)));
		x2(j)=x(Seg3(poselem(j),1))+weight3*(x(Seg3(poselem(j),2))-x(Seg3(poselem(j),1)));
		y1(j)=y(Seg2(poselem(j),1))+weight2*(y(Seg2(poselem(j),2))-y(Seg2(poselem(j),1)));
		y2(j)=y(Seg3(poselem(j),1))+weight3*(y(Seg3(poselem(j),2))-y(Seg3(poselem(j),1)));
		edge_l(j,1)=Seg2_num(poselem(j));
		edge_l(j,2)=Seg3_num(poselem(j));
	else
		%it shoud not go here
	end
end

%now that we have the segments, we must try to connect them...

%loop over the subcontours
contours=struct([]);

while ~isempty(edge_l),

	%take the right edge of the second segment and connect it to the next segments if any
	e1=edge_l(1,1);   e2=edge_l(1,2);
	xc=[x1(1);x2(1)]; yc=[y1(1);y2(1)];

	%erase the lines corresponding to this edge
	edge_l(1,:)=[];
	x1(1)=[]; x2(1)=[];
	y1(1)=[]; y2(1)=[];

	[ro1,co1]=find(edge_l==e1);

	while ~isempty(ro1)

		if co1==1,
			xc=[x2(ro1);xc]; yc=[y2(ro1);yc];

			%next edge:
			e1=edge_l(ro1,2);

		else
			xc=[x1(ro1);xc]; yc=[y1(ro1);yc];

			%next edge:
			e1=edge_l(ro1,1);
		end

		%erase the lines of this
		edge_l(ro1,:)=[];
		x1(ro1)=[]; x2(ro1)=[];
		y1(ro1)=[]; y2(ro1)=[];

		%next connection
		[ro1,co1]=find(edge_l==e1);
	end

	%same thing the other way (to the right)
	[ro2,co2]=find(edge_l==e2);

	while ~isempty(ro2)

		if co2==1,
			xc=[xc;x2(ro2)]; yc=[yc;y2(ro2)];

			%next edge:
			e2=edge_l(ro2,2);
		else
			xc=[xc;x1(ro2)]; yc=[yc;y1(ro2)];

			%next edge:
			e2=edge_l(ro2,1);
		end

		%erase the lines of this
		edge_l(ro2,:)=[];
		x1(ro2)=[]; x2(ro2)=[];
		y1(ro2)=[]; y2(ro2)=[];

		%next connection
		[ro2,co2]=find(edge_l==e2);
	end

	%save xc,yc contour: 
	contours(end+1).x=xc;
	contours(end).y=yc;
	contours(end).name='';
	contours(end).nods=length(xc);
	contours(end).density=1;
	contours(end).closed=0;

end
