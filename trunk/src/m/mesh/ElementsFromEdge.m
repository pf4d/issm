function edgeelements=ElementsFromEdge(elements,A,B)
%ELEMENTSFROMEDGE: find elements connected to one edge defined by nodes A and B
%
%   Usage: edgeelements=ElementsFromEdge(elements,A,B) 
%
%   Eg:    edgeelements=ElementsFromEdge(md.mesh.elements,tip1,tip2)
%
%
edgeelements=find(...
	(elements(:,1)==A & elements(:,2)==B )| ...
	(elements(:,1)==A & elements(:,3)==B )| ...
	(elements(:,2)==A & elements(:,3)==B )| ...
	(elements(:,2)==A & elements(:,1)==B )| ...
	(elements(:,3)==A & elements(:,1)==B )| ...
	(elements(:,3)==A & elements(:,2)==B ));