%TheoryNLOS2.m


clear 
close all;
C=0.007;
mu=2;
Rt=40;
% syms r;
beta=0.01;
beta0=0.1;
omega= pi/2;
p=omega/2/pi;
lambda=2;
R=100;
lamT = 100e-6;
bt= @(r) 2/(C/mu)^2./(Rt^2-r.^2).*(C/mu*(Rt-r)-log((1+C/mu*Rt)./(1+C/mu*r)));


PbNLOS = @(r) (exp(-bt(r)*lambda)-bt(r)*exp(-lambda)).*(r<Rt)+1-(r<Rt);
PbLOS  = @(r) p*exp(-(beta*r+beta0))./(1+C/mu*r);
% temp=@(r) 1-exp(-(beta*r+beta0))./(1+C*r) *2.*r/R^2;
atInt=@(r) PbLOS(r).*PbNLOS(r)*2.*r/R^2;
at=1-integral(atInt,1,100);
pB = exp(-at*lamT*pi*R^2)

iB=1;
iD=1;
c=0.01;
aPrimeExp = @(r) exp(-(beta(iD)*r+beta0(iD)))./(1+c*r) *2.*r/R^2;
aPrime(iB,iD) = integral(aPrimeExp, 0, R);
