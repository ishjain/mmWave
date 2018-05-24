% May23: generate hexagonal cell theoretical results.
%May 20: hexagon
% numerical analysis with hexagonal cells
% Run  sbatch --array=1-6 mybatch.sbatch

clear
close all

aID = getenv('SLURM_ARRAY_TASK_ID')

densityBL = [0.01,0.1,0.2];
densityAP = [50,100,200,300,400,500]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, pi/3, pi/2];

indT = str2num(aID);
lamT = densityAP(indT); %lambda
D = sqrt(2./(3*sqrt(3)*lamT));
ri=[repmat(sqrt(3)*D,1,6) repmat(3*D,1,6) repmat(2*sqrt(3)*D,1,6)];
phi = [pi/6, -pi/6, pi/2, -pi/2, 5*pi/6, -5*pi/6,...
    0, pi/3, -pi/3, 2*pi/3, -2*pi/3, pi,...
    pi/6, -pi/6, pi/2, -pi/2, 5*pi/6, -5*pi/6 ];
di = @(r,theta) sqrt(ri.^2+r.^2-2*ri.*r.*cos(phi-theta));

V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 5;
mu = 2;
frac = (hb-hr)/(ht-hr);
temp = 2/pi*V*frac;
% lamB = 0.1;


% for indT = 1:length(densityAP)
for indB = 1:length(densityBL)
    %         for indO = 1:length(omegaVal)
    
    lamB = densityBL(indB);
    %             omega = omegaVal(indO);
    
    C = temp*lamB;
    
    
    condProb =  @(r,theta) exp(sum(log((C/mu.*di(r,theta))./(1+C/mu.*di(r,theta))).*(di(r,theta)<=100))).*2.*r/D^2 * 1/(2*pi);
    
    results(indB) = integral(@(r)integral(@(theta)condProb(r,theta),0,2*pi,'ArrayValued',true),0,D,'ArrayValued',true);
    
    
    %         end
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
