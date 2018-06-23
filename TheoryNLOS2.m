%TheoryNLOS2.m
%June 22: NLOS model in JSAC

clear
close all;


wannaplot=1;
V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 5;
frac = (hb-hr)/(ht-hr); %fraction depends on heights
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
densityBS = [50,100,200,300,400,500]*10^(-6); %BS=Base Station
densityBL = [.01  0.1]; %Dynamic BLockers
densityD = [1e-9,0.0001]; %D = static blockage
omegaVal = [0 pi/3];

nBS = length(densityBS);
nBL = length(densityBL);
nD = length(densityD);
nO = length(omegaVal);
El = 10; %m
Ew = 10; %m

% !!!-------------------------------!!!
% pB = zeros(nBS,nBL,nO);
% pBCond=zeros(nBS,nBL,nO);
% freq=zeros(nBS,nBL,nO);
% freqCond=zeros(nBS,nBL,nO);
% durCond=zeros(nBS,nBL,nO);



%%NLOS parameters----------
Rt=40;
kappa = 3;

% C=0.007;
% mu=2;

% % syms r;
% beta=0.01;
% beta0=0.1;
% omega= pi/2;
% p=omega/2/pi;
% lambda=2;
% R=100;
% lamT = 100e-6;

for iT = 1:nBS
    tempind = 0;
    for iB = 1:nBL
        for iD = 1:nD
            for iO = 1:nO
                tempind=tempind+1; %increment from zero after every BSdensity value
                lamT = densityBS(iT); %lambda BS
                lamB = densityBL(iB);
                lamD = densityD(iD);
                omega = omegaVal(iO);
                p(iO) =1- omega/(2*pi);
                C(iB) = 2/pi*lamB*V*frac;
                beta(iD) = 2/pi*lamD*(El+Ew);
                beta0(iD) = lamD*El*Ew;
                
                %This required for duration %Please double check!!!!!!!
                q(iD) = 2*exp(-beta0(iD))/(beta(iD)^2*R^2) * ...
                    (1-(1+beta(iD)*R)*exp(-beta(iD)*R));
                %Prepare for coverage
                %                 qt(iD) = 2*exp(-beta0(iD))/(beta(iD)^2*Rt^2) *...
                %                     ((1+beta(iD)*Rt)*exp(-beta(iD)*Rt)-(1+beta(iD)*R)*exp(-beta(iD)*R));
                
                qt(iD) = Rt^2/R^2 - 2*p(iO)*exp(-beta0(iD))/(R^2*beta(iD)^2)*...
                    (exp(-beta(iD)*R)*(1+beta(iD)*R)-exp(-beta(iD)*Rt)*(1+beta(iD)*Rt));
                
                %Define Coverage Probability
                pNoCoverage(iT,iD,iO) = exp(-qt(iD)*lamT*pi*R^2);
                pCoverage(iT,iD,iO) = 1-pNoCoverage(iT,iD,iO);
                pc(iT,tempind) =pCoverage(iT,iD,iO); %temp to check its correcteness
                %Prepare for blockage Probability: numerical integration
                %bt is related to NLOS blockage probability PbNLOS
                bt= @(r) 2/(C(iB)/mu)^2./(Rt^2-r.^2).*(C(iB)/mu*(Rt-r)-log((1+C(iB)/mu*Rt)./(1+C(iB)/mu*r)));
                PbNLOS = @(r) ((exp(-bt(r)*kappa)-bt(r)*exp(-kappa)).*(r<Rt)+1-(r<Rt));
                PbLOS  = @(r) (1-p(iO)*exp(-(beta(iD)*r+beta0(iD)))./(1+C(iB)/mu*r));
                % temp=@(r) 1-exp(-(beta*r+beta0))./(1+C*r) *2.*r/R^2;
                atInt=@(r)  (PbNLOS(r)).*(PbLOS(r))*2.*r/R^2;
                at(iB)=1-integral(atInt,Rt,R,'absTol',1e-100);
                
                %
                
                
                %1. Marginal prob of Blockage (open area)
                pB(iT,tempind) = exp(-at(iB)*lamT*pi*R^2);
                
                
                %2. Conditional prob of blockage given coverage(newly defined)
                pBCond(iT,tempind) = (pB(iT,tempind)-...
                    pNoCoverage(iT,iD,iO))/pCoverage(iT,iD,iO);
                %5. Conditional expectation of duration of bl given coverage
                dur(iT,tempind) = 1/(p(iO)*q(iD)*lamT*pi*R^2+kappa*lamT*pi*Rt^2);
                durCond(iT,tempind) = dur(iT,tempind)*1000/pCoverage(iT,iD,iO);
                
                %%put column title for saving in csv files
                colTitle{1}='lamT';
                if(lamD<1e-7),lamD=0;end
                colTitle{tempind+1} = strcat('lamB',num2str(lamB),...
                    'lamD',num2str(lamD*1e4),'omega',num2str(omega*360/2/pi));
                
                %Put legends
                legendArray{tempind} = strcat(' \lambda_B=',num2str(lamB),...
                    '\lambda_D=',num2str(lamD*1e4),' \omega=',num2str(omega*360/2/pi));
            end
        end
    end
    
end


writetable(cell2table([colTitle; num2cell([densityBS'*10^4, pBCond])]),...
    'figures2/theory_pB_NLOS.csv','writevariablenames',0);
writetable(cell2table([colTitle; num2cell([densityBS'*10^4,durCond])]),...
    'figures2/theory_durCond_NLOS.csv','writevariablenames',0);

if(wannaplot)
    figure(1);grid on;
    semilogy(densityBS,pB);
    ylim([1e-6,1]);
    title('Marginal prob of Blockage')
    legend(legendArray);
    
    
    figure(2);grid on;
    semilogy(densityBS,pBCond);
    title('Conditional prob of Bl given coverage');
    ylim([1e-6,1])
    legend(legendArray);
    
    figure(5); grid on;
    plot(densityBS,durCond)
    title('Conditional expectated duration of bl given coverage')
    legend(legendArray);
    
end

%% coverage plots
clear
wannaplot=1;

R = 100; %m Radius

omegaVal = [0 pi/3];
densityBS = (50:500)*10^-6;
% densityBL = [.01  0.1]; %Dynamic BLockers
densityD = [1e-9,0.0001]; %D = static blockage
omegaVal = [0, pi/3 ];

nBS = length(densityBS);

nO = length(omegaVal);
El = 10; %m
Ew = 10; %m

Rt=40;
kappa = 3;


nD=length(densityD);


for iT=1:nBS
    tempind = 0;
    for iD = 1:nD
        for iO = 1:nO
            
            tempind=tempind+1; %increment from zero after every BSdensity value
            
            lamT = densityBS(iT);
            lamD = densityD(iD);
            omega = omegaVal(iO);
            p(iO) =1- omega/(2*pi);
            %         C(iB) = 2/pi*lamB*V*frac;
            beta(iD) = 2/pi*lamD*(El+Ew);
            beta0(iD) = lamD*El*Ew;
            
            
            qt(iD) = Rt^2/R^2 - 2*p(iO)*exp(-beta0(iD))/(R^2*beta(iD)^2)*...
                (exp(-beta(iD)*R)*(1+beta(iD)*R)-exp(-beta(iD)*Rt)*(1+beta(iD)*Rt));
            
            %Define Coverage Probability
            pNoCoverage(iD,iO) = exp(-qt(iD)*lamT*pi*R^2);
            pCoverage(iT,tempind) = 1-pNoCoverage(iD,iO);
            %%put column title for saving in csv files
            colTitle2{1}='lamT';
            if(lamD<1e-7),lamD=0;end
            colTitle2{tempind+1} = strcat(...
                'lamD',num2str(lamD*1e4),'omega',num2str(omega*360/2/pi));
            legendArray2{tempind} = strcat(...
                '\lambda_D=',num2str(lamD*1e4),' \omega=',num2str(omega*360/2/pi));
        end
    end
end
writetable(cell2table([colTitle2; num2cell([densityBS'*10^4,pCoverage])]),...
    'figures2/coverage_NLOS.csv','writevariablenames',0);

figure(3);grid on;
semilogy(densityBS',pCoverage);
title('coverage NLOS');
% ylim([1e-6,1])
legend(legendArray2);
