function rigidity=buddjacka(temperature)
% BUDDJACKA - calculates ice rigidity as a function of temperature
%
%   rigidity (in s^(1/3)Pa) is the flow law parameter in the flow law sigma=B*e(1/3)
%   Budd and Jacka (1989)
%   temperature is in Kelvin degrees
%
%   Usage:
%      rigidity=buddjacka(temperature)

if any(temperature<0)
	error('input temperature should be in Kelvin (positive)');
end
T=temperature-273.15;

%BJtable=[-5.0000000e-02   1.1690833e-05
%-1.0000000e+00   6.6642879e-06
%-2.0000000e+00   4.0324422e-06
%-5.0000000e+00   1.3461342e-06
%-1.0000000e+01   4.6586675e-07
%-1.5000000e+01   2.2686290e-07
%-2.0000000e+01   1.1855922e-07
%-2.5000000e+01   6.3886499e-08
%-3.0000000e+01   3.5479579e-08
%-3.5000000e+01   1.9228991e-08
%-4.0000000e+01   9.1625910e-09
%-5.0000000e+01   1.4769247e-09];
%
%Temp=BJtable(:,1);
%Ao=BJtable(:,2)*1e-18; %conversion from MPa^-3 to Pa^-3
%Ae=Ao*(2/3)^((1-3)/2);
%B=Ae.^(-1/3);
%fittedmodel=fit(Temp,B,'cubicspline');
%rigidity=fittedmodel(T);
%return

%Temp=[-50:1:0];                                    
%Bcall=buddjacka(Temp+273.15);
%Bbjall=buddjacka(Temp+273.15);
%Acall=(Bcall.^-3)*(2/3)*1e18;
%Abjall=(Bbjall.^-3)*(2/3)*1e18
%semilogy(Temp,Acall,'--k',Temp,Abjall,'k')
%semilogy(Temp,Bcall,'--k',Temp,Bbjall,'k')

rigidity=zeros(size(T));
pos=find(T<=-40);
rigidity(pos)=1e8*(-0.000237326134296*(T(pos)+50).^3+ 0.017054655749852*(T(pos)+50).^2-0.496435984007500*(T(pos)+50)+7.670967258739796);
pos=find(-40<T & T<=-35);
rigidity(pos)=1e8*(-0.000237326134296*(T(pos)+40).^3+ 0.009934871720961*(T(pos)+40).^2-0.226540709299368*(T(pos)+40)+4.174746859353635);
pos=find(-35<T & T<=-30);
rigidity(pos)=1e8*(-0.000293001369376*(T(pos)+35).^3+ 0.006374979706516*(T(pos)+35).^2-0.144991452161983*(T(pos)+35)+3.260749339093782);
pos=find(-30<T & T<=-25);
rigidity(pos)=1e8*(-0.000053702836500*(T(pos)+30).^3+ 0.001979959165871*(T(pos)+30).^2-0.103216757800049*(T(pos)+30)+2.658541399774723);
pos=find(-25<T & T<=-20);
rigidity(pos)=1e8*( 0.000006906867543*(T(pos)+25).^3+ 0.001174416618375*(T(pos)+25).^2-0.087444878878821*(T(pos)+25)+2.185243735358781);
pos=find(-20<T & T<=-15);
rigidity(pos)=1e8*(-0.000015460250554*(T(pos)+20).^3+ 0.001278019631513*(T(pos)+20).^2-0.075182697629382*(T(pos)+20)+1.778243114866868);
pos=find(-15<T & T<=-10);
rigidity(pos)=1e8*(-0.000110386100241*(T(pos)+15).^3+ 0.001046115873209*(T(pos)+15).^2-0.063562020105772*(T(pos)+15)+1.432347586188582);
pos=find(-10<T & T<=-5);
rigidity(pos)=1e8*(-0.000108595885218*(T(pos)+10).^3+-0.000609675630408*(T(pos)+10).^2-0.061379818891767*(T(pos)+10)+1.126892119959808);
pos=find(-5<T & T<=-2);
rigidity(pos)=1e8*( 0.000173187986430*(T(pos)+5).^3+-0.002238613908676*(T(pos)+5).^2-0.075621266587187*(T(pos)+5)+0.791176649088537);
pos=find(-2<T & T<=-1);
rigidity(pos)=1e8*( 0.000429499435151*(T(pos)+2).^3+-0.000679922030808*(T(pos)+2).^2-0.084376874405640*(T(pos)+2)+0.548841399782495);
pos=find(-1<T);
rigidity(pos)=1e8*( 0.000429499435151*(T(pos)+1).^3+ 0.000608576274646*(T(pos)+1).^2-0.084448220161802*(T(pos)+1)+0.464214102781199);

%Now make sure that rigidity is positive
pos=find(rigidity<0);        rigidity(pos)=10^6;
