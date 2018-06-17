%TheoryNLOS
%June 4: Numerical integration 6 times and infinite sum 2 times =D

clear
close all
aID = getenv('SLURM_ARRAY_TASK_ID')
if(isempty(aID))
    warning('aID is empty. Replacing it with 1.')
    aID = '1'; %Runs only for first value of AP density when aID=1
end
aIDnum = str2num(aID);


%Run: sbatch --array=1-24 mybatch.sbatch

%Range of values
densityAP = [50,100,200,300,400,500]*10^(-6);%(1:1:10)/10^4;
densityBL = [0.01,0.1];
densityRef = [10,50]*10^(-6); %Reflector density

for ap = 1:length(densityAP)
    for bl = 1:length(densityBL)
        for ref = 1:length(densityRef)
            %             mylist{ap,bl,ref} = [densityAP(ap),densityBL(bl),densityRef(ref)];
            mylist{ap,bl,ref} = [ap,bl,ref];
        end
    end
end
indT= mylist{aIDnum}(1); %Index of AP
indB= mylist{aIDnum}(2); %Index of Blocker
indL= mylist{aIDnum}(3); %Index of Reflector
%Parameters
V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 5;
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
temp = 2/pi*V*frac*R;

R_NLOS = 60; %m
dmin = 10; %m
dmax = 50; %m





%Arguments
%BS location (r,theta): gives [p1,q1]
%Reflector midpoint locaton (l,phi): gives [p2,q2]
%Reflector length and orientation (d,psi): gives following
%Reflector end points A:[x1,y1], B:[x2,y2]
%User's reflection point O':[a,b]

%Probability distributions
f_r = @(r) 2*r/R^2; %range from [0,R]
f_theta = 1/2/pi; %Unif[0,2pi]
f_l = @(l) 2*l/R^2; %range from [0,R]
f_phi =  1/2/pi; %Unif[0,2pi]
f_d = 1/(dmax-dmin); %Unif[dmin,dmax]
f_psi = 1/2/pi; %Unif[0,2pi]

%End point coordinates:
p1 = @(r,theta) r.*cos(theta);
q1 = @(r,theta) r.*sin(theta);
p2 = @(l,phi) l.*cos(phi);
q2 = @(l,phi) l.*sin(phi);
a = @(l,phi,psi) 2*l.*sin(psi-phi).*sin(psi);
b = @(l,phi,psi) (-2)*l.*sin(psi-phi).*cos(psi);
OA = @(l,phi,d,psi) sqrt((d/2+l.*cos(psi-phi)).^2+(l.*sin(psi-phi)).^2);
OB = @(l,phi,d,psi) sqrt((d/2-l.*cos(psi-phi)).^2+(l.*sin(psi-phi)).^2);
angleOAS = @(l,phi,d,psi) psi + atan((d/2+l.*cos(psi-phi))/(l.*sin(psi-phi)));
angleOBS = @(l,phi,d,psi) psi - atan((d/2-l.*cos(psi-phi))/(l.*sin(psi-phi)));
x1 = @(l,phi,d,psi) OA(l,phi,d,psi).*cos(angleOAS(l,phi,d,psi));
y1 = @(l,phi,d,psi) OA(l,phi,d,psi).*sin(angleOAS(l,phi,d,psi));
x2 = @(l,phi,d,psi) OB(l,phi,d,psi).*cos(angleOBS(l,phi,d,psi));
y2 = @(l,phi,d,psi) OB(l,phi,d,psi).*sin(angleOBS(l,phi,d,psi));

%NLOS coverage indicator random variable
%We need to add self blockage???

d_NLOS = @(r,theta,l,phi,psi) sqrt(r.^2+4*l.^2*sin(psi-phi).^2+4*r.*l.*sin(psi-phi)*sin(theta-psi));
d_NLOSapprox = @(r,l) r+2*l;
% r=50;
% d=10;
% theta=pi/4;

%Indicates whether NLOS path exists (1) or not (0)
IndicatorRV= @(r,theta,l,phi,d,psi) ((q1(r,theta)-b(l,phi,psi)-...
    ((y1(l,phi,d,psi)-b(l,phi,psi))*(p1(r,theta)-a(l,phi,psi)))/(x1(l,phi,d,psi)-a(l,phi,psi))) >=0) && ...
    ((q2(l,phi)-b(l,phi,psi)-...
    ((y2(l,phi,d,psi)-b(l,phi,psi))*(p2(l,phi)-a(l,phi,psi)))/(x2(l,phi,d,psi)-a(l,phi,psi))) >=0) &&...
    (d_NLOS(r,theta,l,phi,psi)<=R_NLOS);

% for indB = 1:length(densityBL)
%avoid for loop :) using 24 parallel computation
%     lamB = densityBL(indB);


lamB = densityBL(indB);
C(indB) = 2/pi*lamB*V*frac;

%Conditional NLOS Probability conditioned on (r,theta,l,phi,d,psi,nT,nB)
P_NLOS_Cond = @(r,theta,l,phi,d,psi)(C(indB)/mu*d_NLOSapprox(r,l)/...
    (1+C(indB)/mu*d_NLOSapprox(r,l)))*...
    IndicatorRV(r,theta,l,phi,d,psi) +1-IndicatorRV(r,theta,l,phi,d,psi);


int1 = @(r,theta,l,phi,d) integral(@(psi)P_NLOS_Cond(r,theta,l,phi,d,psi)*f_psi,0,2*pi,'ArrayValued',true);
int2 = @(r,theta,l,phi) integral(@(d)int1(r,theta,l,phi,d).*f_d,dmin,dmax,'ArrayValued',true);
int3 = @(r,theta,l) integral(@(phi)int2(r,theta,l,phi).*f_phi,0,2*pi,'ArrayValued',true);
int4 = @(r,theta) integral(@(l)int3(r,theta,l).*f_l(l),0,R,'ArrayValued',true);
int5 = @(r) integral(@(theta)int4(r,theta)*f_theta,0,2*pi,'ArrayValued',true);

tic 
a_NLOS(indB) = integral(@(r)int5(r).*f_r(r),0,R,'ArrayValued',true);

toc
tic

%Final P_NLOS is obtained as a function of lamT,lamB,lamL
%     for indT = 1:length(densityAP)
%         for indL = 1:length(densityRef)
lamT = densityAP(indT); %lambda
lamL = densityRef(indL); %lambda RefLector
syms nL nT
P_NLOS_temp = symsum(  a_NLOS(indB).^(nL.*nT).*exp(-lamL).*lamL.^(nL)./factorial(nL), nL, 0, Inf);
P_NLOS = vpa(symsum( P_NLOS_temp.*exp(-lamT).*lamT.^(nT)./factorial(nT), nT, 0, Inf));
toc
%         end
%     end
% end
% %

csvwrite(strcat('output',num2str(aID),'.csv'),P_NLOS)
