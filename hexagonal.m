% May23: generate hexagonal cell theoretical results.
%May 20: hexagon
% numerical analysis with hexagonal cells
% HowtoRunOnHPC?:  sbatch --array=1-6 mybatch.sbatch %for 6 values of AP density

clear
close all

aID = getenv('SLURM_ARRAY_TASK_ID')
if(isempty(aID))
  warning('aID is empty. Replacing it with 1.')  
  aID = '1'; %Runs only for first value of AP density when aID=1
end
densityBL = [0.01,0.1];
densityAP = [50,100,200,300,400,500]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, pi/3];

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
hb = 1.8;
hr = 1.4;
ht = 5;
mu = 2;
R=100;
frac = (hb-hr)/(ht-hr);
temp = 2/pi*V*frac;
% lamB = 0.1;


% for indT = 1:length(densityAP)
%iterate over BS density
for indO = 1:length(omegaVal)
    if(indO==1)
        di = @(r,theta) sqrt(ri.^2+r.^2-2*ri.*r.*cos(phi-theta));
    elseif(indO==2)
        di = @(r,theta) sqrt(ri_worst60.^2+r.^2-2*ri_worst60.*r.*cos(phi_worst60-theta));
    end
    %iterate over blocker density
    for indB = 1:length(densityBL)
        lamB = densityBL(indB);
        C = temp*lamB;
        condProb =  @(r,theta) exp(sum(log((C/mu.*di(r,theta))./(1+C/mu.*di(r,theta))).*(di(r,theta)<=R))).*...
            2.*r/d^2 * 1/(2*pi);
        tic
        results(indB,indO) = integral(@(r)integral(@(theta)condProb(r,theta),0,2*pi,'ArrayValued',true),...
            0,d,'ArrayValued',true);        
        toc
    end
    
end
% end
csvwrite(strcat('output',num2str(aID),'.csv'),results)
%
% % condProb = @(r,theta,D) r .*cos(theta)
%
% condProb =  @(r,theta) exp(sum(log((C/mu.*di(r,theta))./(1+C/mu.*di(r,theta))).*(di(r,theta)<=100))).*2.*r/D^2 * 1/(2*pi);
%
% results = integral(@(r)integral(@(theta)condProb(r,theta),0,2*pi,'ArrayValued',true),0,D,'ArrayValued',true)
%
