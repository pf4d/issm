function [x,y] = ll2xy(lat,lon,sgn,central_meridian,standard_parallel)  
%LL2XY - converts lat long to polar stereographic
%
%   Converts from geodetic latitude and longitude to Polar 
%   Stereographic (X,Y) coordinates for the polar regions.
%   Author: Michael P. Schodlok, December 2003 (map2ll)
%
%   Usage:
%      [x,y] = ll2xy(lat,lon,sgn)
%      [x,y] = ll2xy(lat,lon,sgn,central_meridian,standard_parallel)
%
%      - sgn = Sign of latitude +1 : north latitude (default is mer=45 lat=70)
%                               -1 : south latitude (default is mer=0  lat=71)

%Get central_meridian and standard_parallel depending on hemisphere
if nargin==5,
	delta = central_meridian;
	slat  = standard_parallel;
elseif nargin==3
	if sgn == 1,
		delta = 45; slat = 70;
		disp('Info: creating coordinates in polar stereographic (Std Latitude: 70ºN Meridian: 45º)');
	elseif sgn==-1,
		delta = 0;  slat = 71;
		disp('Info: creating coordinates in polar stereographic (Std Latitude: 71ºS Meridian: 0º)');
	else
		error('Sign should be either +1 or -1');
	end
else
	help ll2xy
	error('bad usage');
end

% Conversion constant from degrees to radians
cde  = 57.29577951;
% Radius of the earth in meters
re   = 6378.273*10^3;
% Eccentricity of the Hughes ellipsoid squared
ex2   = .006693883;
% Eccentricity of the Hughes ellipsoid
ex    =  sqrt(ex2);

latitude  = abs(lat) * pi/180.;
longitude = (lon + delta) * pi/180.;

% compute X and Y in grid coordinates.
T = tan(pi/4-latitude/2) ./ ((1-ex*sin(latitude))./(1+ex*sin(latitude))).^(ex/2);

if (90 - slat) <  1.e-5 
	rho = 2.*re*T/sqrt((1.+ex)^(1.+ex)*(1.-ex)^(1.-ex));
else
	sl  = slat*pi/180.;
	tc  = tan(pi/4.-sl/2.)/((1.-ex*sin(sl))/(1.+ex*sin(sl)))^(ex/2.);
	mc  = cos(sl)/sqrt(1.0-ex2*(sin(sl)^2));
	rho = re*mc*T/tc;
end

y = -rho .* sgn .* cos(sgn.*longitude);
x =  rho .* sgn .* sin(sgn.*longitude);

[cnt1,cnt2] = find(latitude(:) >= pi / 2.);

if cnt1
	x(cnt1) = 0.0;
	y(cnt1) = 0.0;
end
