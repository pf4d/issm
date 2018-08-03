function expswapxy(filename)
%EXPSWAP - swap x and y fields
% 
%   Usage:
%      expswap(file)
%
%   See also EXPMASTER, EXPDOC

contours=expread(filename,1);

newcontours=contours(1);

for i=1:length(contours), 
	contour=contours(i);
	newcontour=contour;
	newcontour.x=contour.y;
	newcontour.y=contour.x;
	newcontours(end+1)=newcontour;
end
newcontours=newcontours(2:end);

expwrite(newcontours,filename);
