function plot_elementstype(md,options,width,i)
%PLOT_ELEMENTSTYPE - plot elements type
%
%   Usage:
%      plot_elementstype(md,options,width,i);
%
%   See also: PLOTMODEL

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
[data datatype]=processdata(md,md.flowequation.element_equation,options);

%edgecolor?
edgecolor=getfieldvalue(options,'edgecolor','none');

%plot
subplot(width,width,i);

if is2d
	pos=find(data==0);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	p1=patch( 'Faces', [A B C], 'Vertices', [x y z],'CData',0,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==1);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	p2=patch( 'Faces', [A B C], 'Vertices', [x y z],'CData',2,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==2);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	p3=patch( 'Faces', [A B C], 'Vertices', [x y z],'CData',3,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==3);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	p4=patch( 'Faces', [A B C], 'Vertices', [x y z],'CData',4,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==4);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	p5=patch( 'Faces', [A B C], 'Vertices', [x y z],'CData',5,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==5);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	p6=patch( 'Faces', [A B C], 'Vertices', [x y z],'CData',6,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==6);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	p7=patch( 'Faces', [A B C], 'Vertices', [x y z],'CData',7,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==7);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	p8=patch( 'Faces', [A B C], 'Vertices', [x y z],'CData',8,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==8);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	p9=patch( 'Faces', [A B C], 'Vertices', [x y z],'CData',9,'FaceColor','flat','EdgeColor',edgecolor);
else
	pos=find(data==0);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	p1=patch( 'Faces', [A B C],'Vertices', [x y z],'CData', 0,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', 0,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', 0,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', 0,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', 0,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==1);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	p2=patch( 'Faces', [A B C],'Vertices', [x y z],'CData', 1,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', 1,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', 1,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', 1,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', 1,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==2);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	p3=patch( 'Faces', [A B C],'Vertices', [x y z],'CData', 2,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', 2,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', 2,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', 2,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', 2,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==3);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	p4=patch( 'Faces', [A B C],'Vertices', [x y z],'CData', 3,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', 3,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', 3,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', 3,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', 3,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==3);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	p5=patch( 'Faces', [A B C],'Vertices', [x y z],'CData', 3,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', 3,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', 3,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', 3,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', 3,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==4);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	p6=patch( 'Faces', [A B C],'Vertices', [x y z],'CData', 4,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', 4,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', 4,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', 4,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', 4,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==5);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	p7=patch( 'Faces', [A B C],'Vertices', [x y z],'CData', 5,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', 5,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', 5,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', 5,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', 5,'FaceColor','flat','EdgeColor',edgecolor);
	%HOFS elements
	pos=find(data==7);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	p8=patch( 'Faces', [A B C],'Vertices', [x y z],'CData', 7,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', 7,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', 7,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', 7,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', 7,'FaceColor','flat','EdgeColor',edgecolor);
	pos=find(data==6);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	p9=patch( 'Faces', [A B C],'Vertices', [x y z],'CData', 6,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', 6,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', 6,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', 6,'FaceColor','flat','EdgeColor',edgecolor);
	patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', 6,'FaceColor','flat','EdgeColor',edgecolor);
end
legend([p1 p2 p3 p4 p5 p6 p7 p8 p9],...
		'None','SIA','SSA','L1L2','HO',...
		'SSAHO','FS','SSAFS','HOFS');

%apply options
options=addfielddefault(options,'title','Elements type');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
