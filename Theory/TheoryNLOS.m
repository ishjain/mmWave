% Written by Ish Jain
% NYU Tandon School of Engineering
% Date: June 2018
%
% Descripton:
% Generate theoretical plots based on NLOS model presented in JSAC paper


clear
close all;

wannaplotCoverage=0; %change only when focussing on coverage
wannaplotCellRadius=0; %change only when focussing on cell radius
wannaplot=1;
wannaSaveFiles=0; %Saving mat files for plotting later

V = 1; %velocity m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr); %fraction depends on heights
mu = 2; %Expected bloc dur =1/mu
R = 200; %m Radius
densityBS = [50,100,200,300,400,500]*10^(-6); %BS=Base Station
densityBL = [0.01, 0.1  ]; %Dynamic BLockers
densityD = [1e-9,0.0001]; %D = static blockage
omegaVal = [0, pi/3];

nBS = length(densityBS);
nBL = length(densityBL);
nD = length(densityD);
nO = length(omegaVal);

%Expected length and width of buildings
El = 10; %m
Ew = 10; %m


%%NLOS parameters----------

gammaNLOS = 5;
PLE = 2.69;
Rt=R*10^(-gammaNLOS/(10*PLE)); %65; %R^(0.91); %0.91=PLE(LOS)/PLE(NLOS), where PLE=path loss exponent
kappa = 2;


for iT = 1:nBS
    tempind = 0;
    for iB = 1:nBL
        for iD = 1:nD
            for iO = 1:nO
                tempind=tempind+1; %increment from zero after every BSdensity value
                lamT = densityBS(iT); %lambda BS
                lamB = densityBL(iB);
                lamD = densityD(iD);
                omega = omegaVal(iO); %self-blockage angle
                p(iO) =1- omega/(2*pi); %probability of no self-blockage
                C(iB) = 2/pi*lamB*V*frac; %C is defined in paper
                beta(iD) = 2/pi*lamD*(El+Ew); %static blockage parameter
                beta0(iD) = lamD*El*Ew; %static blockage parameter
                
                %This required for duration %Please double check!!!!!!!
                q(iD) = 2*exp(-beta0(iD))/(beta(iD)^2*R^2) * ...
                    (1-(1+beta(iD)*R)*exp(-beta(iD)*R));
                %Prepare for coverage
                %                 qt(iD) = 2*exp(-beta0(iD))/(beta(iD)^2*Rt^2) *...
                %                     ((1+beta(iD)*Rt)*exp(-beta(iD)*Rt)-(1+beta(iD)*R)*exp(-beta(iD)*R));
                
                qt(iD,iO) = Rt^2/R^2 - 2*p(iO)*exp(-beta0(iD))/(R^2*beta(iD)^2)*...
                    (exp(-beta(iD)*R)*(1+beta(iD)*R)-exp(-beta(iD)*Rt)*(1+beta(iD)*Rt));
                
                %Define Coverage Probability
                pNoCoverage(iT,iD,iO) = exp(-qt(iD,iO)*lamT*pi*R^2);
                pCoverage(iT,iD,iO) = 1-pNoCoverage(iT,iD,iO);
                pc(iT,tempind) =pCoverage(iT,iD,iO); %temp to check its correcteness
                %Prepare for blockage Probability: numerical integration
                %bt is related to NLOS blockage probability PbNLOS
                bt= @(r) 2/(C(iB)/mu)^2./(Rt^2-r.^2).*...
                    (C(iB)/mu*(Rt-r)-log((1+C(iB)/mu*Rt)./(1+C(iB)/mu*r))).*(r<=Rt);
                PbNLOS = @(r) ((exp(-bt(r)*kappa)-bt(r)*exp(-kappa)));
                
                
                
                
                %Get LOS blockage probility
                PbLOS  = @(r) (1-p(iO)*exp(-(beta(iD)*r+beta0(iD)))./(1+C(iB)/mu*r));
                % temp=@(r) 1-exp(-(beta*r+beta0))./(1+C*r) *2.*r/R^2;
                
                %Combine LOS and NLOS blockage probability
                atInt=@(r)  (PbNLOS(r)).*(PbLOS(r))*2.*r/R^2;
                at(iB)=1-integral(atInt,0,Rt-1e-20,'absTol',1e-100)-integral(atInt,Rt+1e-20,R,'absTol',1e-100);
                
                
                %1. Marginal prob of Blockage (open area)
                pB(iT,tempind) = exp(-at(iB)*lamT*pi*R^2);
                
                
                %2. Conditional prob of blockage given coverage(newly defined)
                pBCond(iT,tempind) = (pB(iT,tempind)-...
                    pNoCoverage(iT,iD,iO))/pCoverage(iT,iD,iO);
                
                %Let's get a lower bound on NLOS blockage probability
                btMin = @(r) 1./(1+C(iB)/mu).*(r<Rt);
                PbNLOSmin = @(r) ((exp(-btMin(r)*kappa)-btMin(r)*exp(-kappa)));
                atIntMin = @(r)  (PbNLOSmin(r)).*(PbLOS(r))*2.*r/R^2;
                atMin(iB) = 1-integral(atIntMin,0,Rt-1e-20,'absTol',1e-100)-integral(atIntMin,Rt+1e-20,R,'absTol',1e-100);
                pBmin(iT,tempind) = exp(-atMin(iB)*lamT*pi*R^2);
                pBCondMin(iT,tempind) = (pBmin(iT,tempind)-...
                    pNoCoverage(iT,iD,iO))/pCoverage(iT,iD,iO);
                
                %5. Conditional expectation of duration of bl given coverage
                dur(iT,tempind) = 1/mu * 1/(p(iO)*q(iD)*lamT*pi*R^2+kappa*lamT*pi*Rt^2);
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

if(wannaSaveFiles)
    writetable(cell2table([colTitle; num2cell([densityBS'*10^4, pBCond])]),...
        'figures/theory_pB_NLOS.csv','writevariablenames',0);
    writetable(cell2table([colTitle; num2cell([densityBS'*10^4, pBCondMin])]),...
        'figures/theory_pB_NLOS_Min.csv','writevariablenames',0);
    writetable(cell2table([colTitle; num2cell([densityBS'*10^4,durCond])]),...
        'figures/theory_durCond_NLOS.csv','writevariablenames',0);
end
if(wannaplot)
    figure(1);grid on;
    semilogy(densityBS,pB);
    
    ylim([1e-6,1]);
    title('Marginal prob of Blockage')
    legend(legendArray);
    
    
    figure(2);grid on;
    semilogy(densityBS,pBCond);
    hold on;
    semilogy(densityBS,pBCondMin);
    title('Conditional prob of Bl given coverage (lower bound also shown)');
    ylim([1e-6,1])
    legend(legendArray);
    %
    figure(5); grid on;
    plot(densityBS,durCond)
    title('Conditional expectated duration of bl given coverage')
    legend(legendArray);
    
end

%% coverage plots
if (wannaplotCoverage)
    clear
    % wannaplot=1;
    
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
    
    Rt=66;
    % kappa = 10;
    
    
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
                
                
                qt(iD,iO) = Rt^2/R^2 - 2*p(iO)*exp(-beta0(iD))/(R^2*beta(iD)^2)*...
                    (exp(-beta(iD)*R)*(1+beta(iD)*R)-exp(-beta(iD)*Rt)*(1+beta(iD)*Rt));
                
                %Define Coverage Probability
                pNoCoverage(iD,iO) = exp(-qt(iD,iO)*lamT*pi*R^2);
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
        'figures/coverage_NLOS.csv','writevariablenames',0);
    
    figure(3);grid on;
    semilogy(densityBS',pCoverage);
    title('coverage NLOS');
    % ylim([1e-6,1])
    legend(legendArray2);
end


%% Effect of cell radius.. same code copied with some modification
if(~exist('wannaplotCellRadius') )
    wannaplotCellRadius=0;
end
if(wannaplotCellRadius)
    
    clear
    % close all
    V = 1; %velocity m/s
    hb = 1.8;
    hr = 1.4;
    ht = 5;
    frac = (hb-hr)/(ht-hr); %fraction depends on heights
    mu = 2; %Expected bloc dur =1/mu
    %     Rvalues = 50:500; %m Radius
    Rvalues = [100, 200];
    
    %     densityBS = [100]*10^(-6); %BS=Base Station
    densityBS = [0.01,1,50,100,150,200,250,300,350,400]*10^(-6); %BS=Base Station
    
    densityBL = [ 0.01  ]; %Dynamic BLockers
    densityD = [0.0001]; %D = static blockage
    omegaVal = [ pi/3];
    
    nR = length(Rvalues);
    nBS = length(densityBS);
    nBL = length(densityBL);
    nD = length(densityD);
    nO = length(omegaVal);
    El = 10; %m
    Ew = 10; %m
    
    %%NLOS parameters----------
 
    PLE = 2.69;
    gammaNLOS = 5;
    Rtvalues = Rvalues*10^(gammaNLOS/(10*PLE));
    kappa = 3;
    for iT = 1:nBS
        tempind = 0;
        for iR = 1:nR
            R = Rvalues(iR);
            Rt= Rtvalues(iR);
            
            
            
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
                        
                        qt(iD,iO) = Rt^2/R^2 - 2*p(iO)*exp(-beta0(iD))/(R^2*beta(iD)^2)*...
                            (exp(-beta(iD)*R)*(1+beta(iD)*R)-exp(-beta(iD)*Rt)*(1+beta(iD)*Rt));
                        
                        %Define Coverage Probability
                        pNoCoverage(iT,iD,iO) = exp(-qt(iD,iO)*lamT*pi*R^2);
                        pCoverage(iT,iD,iO) = 1-pNoCoverage(iT,iD,iO);
                        pc(iT,tempind) =pCoverage(iT,iD,iO); %temp to check its correcteness
                        %Prepare for blockage Probability: numerical integration
                        %bt is related to NLOS blockage probability PbNLOS
                        bt= @(r) 2/(C(iB)/mu)^2./(Rt^2-r.^2).*...
                            (C(iB)/mu*(Rt-r)-log((1+C(iB)/mu*Rt)./(1+C(iB)/mu*r))).*(r<=Rt);
                        PbNLOS = @(r) ((exp(-bt(r)*kappa)-bt(r)*exp(-kappa)));
                        PbLOS  = @(r) (1-p(iO)*exp(-(beta(iD)*r+beta0(iD)))./(1+C(iB)/mu*r));
                        % temp=@(r) 1-exp(-(beta*r+beta0))./(1+C*r) *2.*r/R^2;
                        atInt=@(r)  (PbNLOS(r)).*(PbLOS(r))*2.*r/R^2;
                        at(iB)=1-integral(atInt,0,Rt-1e-6,'absTol',1e-100)-integral(atInt,Rt+1e-6,R,'absTol',1e-100);
                        
                        %1. Marginal prob of Blockage (open area)
                        pB(iT,tempind) = exp(-at(iB)*lamT*pi*R^2);
                        
                        
                        
                        %2. Conditional prob of blockage given coverage(newly defined)
                        pBCond(iR,tempind) = (pB(iT,tempind)-...
                            pNoCoverage(iT,iD,iO))/pCoverage(iT,iD,iO);
                        
                        %Let's get a lower bound on NLOS blockage probability
                        btMin = @(r) 1./(1+C(iB)/mu).*(r<Rt);
                        PbNLOSmin = @(r) ((exp(-btMin(r)*kappa)-btMin(r)*exp(-kappa)));
                        atIntMin = @(r)  (PbNLOSmin(r)).*(PbLOS(r))*2.*r/R^2;
                        atMin(iB) = 1-integral(atIntMin,0,Rt-1e-20,'absTol',1e-100)-integral(atIntMin,Rt+1e-20,R,'absTol',1e-100);
                        pBmin(iT,tempind) = exp(-atMin(iB)*lamT*pi*R^2);
                        pBCondMin(iT,tempind) = (pBmin(iT,tempind)-...
                            pNoCoverage(iT,iD,iO))/pCoverage(iT,iD,iO);
                   
                        %%put column title for saving in csv files
                        colTitle{1}='Radius';
                        if(lamD<1e-7),lamD=0;end
                        colTitle{tempind+1} = strcat('lamB',num2str(lamB),...
                            'lamD',num2str(lamD*1e4),'omega',num2str(omega*360/2/pi),...
                            'R',num2str(R));
                        
                        %Put legends
                        legendArray{tempind} = strcat(' \lambda_B=',num2str(lamB),...
                            '\lambda_D=',num2str(lamD*1e4),' \omega=',num2str(omega*360/2/pi),...
                            'R=',num2str(R));
                    end
                end
            end
        end
    end
    %     writetable(cell2table([colTitle; num2cell([Rvalues', pBCondMin])]),...
    %         'figures/theory_withR_NLOS.csv','writevariablenames',0);
    %     save('figures/R_NLOS.mat','pBCondMin');

    
    figure(8);grid on;
    semilogy(densityBS,pBCondMin, 'LineWidth',2);
    xlabel('BS Density');
    title('LOS Blockage Probability');
    ylim([1e-6,1])
    legend(legendArray);
    %     legend('R=100m','R=200m','R=500m','R=1000m');
end