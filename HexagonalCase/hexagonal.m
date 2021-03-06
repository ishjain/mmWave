% Written by Ish Jain
% NYU Tandon School of Engineering
% Date: June 2018
%
% Description:
% Hexagonal cell case for open park scenario. get theiretical results using
% numerical analysis
%
% Since the code takes lot of time we only run for single BS density in one 
% go or run in parallel for multiple BSs using NYU HPC.
%
% HowtoRunOnHPC?:  sbatch --array=1-12 mybatch.sbatch %for 12 values of AP density

clear
close all


aID = getenv('SLURM_ARRAY_TASK_ID')
if(isempty(aID))
  warning('aID is empty. Replacing it with 1.')  
  aID = '1'; %Runs only for first value of AP density when aID=1
end
densityBL = [0.01];
densityAP = [50,100,150,175,200,225,250,300,350,400,450,500]*10^(-6);%(1:1:10)/10^4;
omegaVal = [ pi/3]; %It will work only for {0,pi/3}, so don't change :P

indT = str2num(aID);
lamT = densityAP(indT); %lambda
d = sqrt(2./(3*sqrt(3)*lamT));
ri=[0,repmat(sqrt(3)*d,1,6) repmat(3*d,1,6) repmat(2*sqrt(3)*d,1,6),...
    repmat(sqrt(21)*d,1,12),repmat(3*sqrt(3)*d,1,6)];

phi = [0, pi/6, -pi/6, pi/2, -pi/2, 5*pi/6, -5*pi/6,...
    0, pi/3, -pi/3, 2*pi/3, -2*pi/3, pi,...
    pi/6, -pi/6, pi/2, -pi/2, 5*pi/6, -5*pi/6,...
    0.06*pi, -0.06*pi, pi/3+0.06*pi, pi/3-0.06*pi, -pi/3+0.06*pi, -pi/3-0.06*pi,...
    2*pi/3+0.06*pi, 2*pi/3-0.06*pi, -2*pi/3+0.06*pi, -2*pi/3-0.06*pi,...
    pi-0.06*pi, -pi+0.06*pi,...
    pi/6, -pi/6, pi/2, -pi/2, 5*pi/6, -5*pi/6];

%The case of worst 60 degree
ri_worst60 = [repmat(sqrt(3)*d,1,4) repmat(3*d,1,5) repmat(2*sqrt(3)*d,1,4),...
    repmat(sqrt(21)*d,1,10),repmat(3*sqrt(3)*d,1,4)];
phi_worst60 = [pi/2, -pi/2, 5*pi/6, -5*pi/6,...
    pi/3, -pi/3, 2*pi/3, -2*pi/3, pi,...
    pi/2, -pi/2, 5*pi/6, -5*pi/6,...
    pi/3+0.06*pi, pi/3-0.06*pi, -pi/3+0.06*pi, -pi/3-0.06*pi,...
    2*pi/3+0.06*pi, 2*pi/3-0.06*pi, -2*pi/3+0.06*pi, -2*pi/3-0.06*pi,...
    pi-0.06*pi, -pi+0.06*pi,...
    pi/2, -pi/2, 5*pi/6, -5*pi/6];

V = 1; %velocity m/s
%height of blocker, UE and BS
hb = 1.8;
hr = 1.4;
ht = 5;
mu = 2;
R=200;
frac = (hb-hr)/(ht-hr);
temp = 2/pi*V*frac;

% iterate over BS density
for indO = 1:length(omegaVal)
    if(indO==1) %omega=0;
        di = @(r,theta) sqrt(ri.^2+r.^2-2*ri.*r.*cos(phi-theta));
    elseif(indO==2) %omega=60degree
        di = @(r,theta) sqrt(ri_worst60.^2+r.^2-2*ri_worst60.*r.*cos(phi_worst60-theta));
    end
    %iterate over blocker density
    for indB = 1:length(densityBL)
        lamB = densityBL(indB);
        C = temp*lamB;
        condProb =  @(r,theta) exp(sum(log((C/mu.*di(r,theta))./(1+C/mu.*di(r,theta))).*(di(r,theta)<=R))).*...
            2.*r/d^2 * 1/(2*pi);
        tic
        results(indB,indO) = integral(@(r)integral(@(theta)condProb(r,theta),0,2*pi,'ArrayValued',true,'AbsTol',1e-50),...
            0,d,'ArrayValued',true,'AbsTol',1e-50);        
        toc
    end
    
end

csvwrite(strcat('output',num2str(aID),'.csv'),results)
