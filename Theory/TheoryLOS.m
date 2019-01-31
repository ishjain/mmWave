% Written by Ish Jain
% NYU Tandon School of Engineering
% Date: June 2018
%
% Description:
%
% Use generalized model and find dynamic blockage probability given
% coverage. Coverage is defined as "atleast one BS out side blocked zone
% (defined by static and self-blockage)"

clear
close all

wannaplotCellRadius=0; %change only when focussing on cell radius
wannaSaveFiles=1; %to save csv files for plotting using Latex.
wannaplot=1; %plot results to visualize

V = 1; %velocity of blocker m/s
hb = 1.8; %height of blocker
hr = 1.4; % height UE
ht = 5;%height BS
frac = (hb-hr)/(ht-hr); %fraction depends on heights
mu = 2; %Expected bloc dur =1/mu

densityBS = [50,100,200,300,400,500]*10^(-6); %BS=Base Station
densityBL = [0.01, 0.1]; %Dynamic BLockers
densityD = [1e-9,0.0001]; %D = static blockage
omegaVal = [0 pi/3]; %self-blockage angle

nBS = length(densityBS);
nBL = length(densityBL);
nD = length(densityD);
nO = length(omegaVal);
El = 10; %building length %m
Ew = 10; %building width %m


pB = zeros(nBS,nBL,nO);
pBCond=zeros(nBS,nBL,nO);

freqCond=zeros(nBS,nBL,nO);
durCond=zeros(nBS,nBL,nO);


R = 200; %m Radius
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
                a(iB) = 2*mu/(R*C(iB))-2*mu^2/(R^2*C(iB)^2)*log(1+R*C(iB)/mu);
                
                % Coverage probability for open park like area
                pnn0(iT,iO)=exp(-p(iO)*lamT*pi*R^2); %prob that n not equal to 0 (pnn0)
                
                % Define parameters for static blockage
                beta(iD) = 2/pi*lamD*(El+Ew);
                beta0(iD) = lamD*El*Ew;
                
                b(iD) = 2*exp(-beta0(iD))/(beta(iD)^2*R^2) * ...
                    (1-(1+beta(iD)*R)*exp(-beta(iD)*R));
                
                %write code to verify b: verified. Both matches
                bExp = @(r) exp(-(beta(iD).*r+beta0(iD))) *2.*r/R^2;
                b(iD) = integral(bExp,0,R);
                
                %coverage probability given by static and self blockage
                pCoverage(iT,iD,iO) = 1-exp(-b(iD)*p(iO)*lamT*pi*R^2);
                pNoCoverage(iT,iD,iO) = 1-pCoverage(iT,iD,iO);
                
                %1. Marginal prob of Blockage (open area)
                pB(iT,iB,iO) = exp(-a(iB)*p(iO)*lamT*pi*R^2);
                
                %1.1 considering static blockage too.
                c=C(iB)/mu;
                
                aPrimeExp = @(r) exp(-(beta(iD)*r+beta0(iD)))./(1+c*r) *2.*r/R^2;
                aPrime(iB,iD) = integral(aPrimeExp, 0, R);
                
                pBprime(iT,tempind) = exp(-aPrime(iB,iD)*p(iO)*lamT*pi*R^2);
                
                
                %2. Conditional prob of blockage given coverage(newly defined)
                pBCondPrime(iT,tempind) = (pBprime(iT,tempind)-...
                    pNoCoverage(iT,iD,iO))/pCoverage(iT,iD,iO);
                
                %3. Expected Freq of blockage
                freqPrime(iT,tempind) = mu*(1-aPrime(iB,iD))*p(iO)*b(iD)...
                    *lamT*pi*R^2*exp(-aPrime(iB,iD)*p(iO)*b(iD)*lamT*pi*R^2);
                
                %4. Conditional expectation of freq of bl given n!=0
                freqCondPrime(iT,tempind) = freqPrime(iT,tempind)/pCoverage(iT,iD,iO);
                
                %5. Conditional expectation of duration of bl given n!=0
                
                %This is used to get duration
                tempb(iT,iD,iO) = p(iO)*b(iD)*lamT*pi*R^2;%*angFrac;
                
                Ei(iT,iB,iO) = ei(tempb(iT,iD,iO))-log(tempb(iT,iD,iO))-0.5772; %exponential integral function
                
                durCondPrime(iT,tempind) = (1-pCoverage(iT,iD,iO))*Ei(iT,iB,iO)/(mu*pCoverage(iT,iD,iO));
                durCondPrime(iT,tempind) = durCondPrime(iT,tempind)*1000;
                %Try different approximations
                durCond(iT,iB,iO) = (1/tempb(iT,iD,iO)+1/tempb(iT,iD,iO)^2)/...
                    (mu*(1-pnn0(iT,iO)));
                durCond2(iT,iB,iO) = (1/tempb(iT,iD,iO))/...
                    (mu*(1-pnn0(iT,iO)));
                durCond2Prime(iT,tempind) = (1/tempb(iT,iD,iO))/...
                    (mu*pCoverage(iT,iD,iO));
                durCond2Prime(iT,tempind)=durCond2Prime(iT,tempind)*1000;
                
                colTitle{1}='lamT';
                if(lamD<1e-7),lamD=0;end
                
                colTitle{tempind+1} = strcat('lamB',num2str(lamB),...
                    'lamD',num2str(lamD*1e4),'omega',num2str(omega*360/2/pi));
                legendArray{tempind} = strcat(' \lambda_B=',num2str(lamB),...
                    '\lambda_D=',num2str(lamD*1e4),' \omega=',num2str(omega*360/2/pi));
                %                 sprintf('\lambda_B=%f,\lambda_D=%d,\omega=%d',lamB,lamD*1e6,omega*360/2/pi);
            end
        end
    end
end


if(wannaSaveFiles)
    writetable(cell2table([colTitle; num2cell([densityBS'*10^4,pBprime])]),...
        'figures/theory_pB.csv','writevariablenames',0)
    writetable(cell2table([colTitle; num2cell([densityBS'*10^4,pBCondPrime])]),...
        'figures/theory_pBCond.csv','writevariablenames',0)
    writetable(cell2table([colTitle; num2cell([densityBS'*10^4,freqPrime])]),...
        'figures/theory_freq.csv','writevariablenames',0)
    writetable(cell2table([colTitle; num2cell([densityBS'*10^4,freqCondPrime])]),...
        'figures/theory_freqCond.csv','writevariablenames',0)
    writetable(cell2table([colTitle; num2cell([densityBS'*10^4,durCond2Prime])]),...
        'figures/theory_durCond.csv','writevariablenames',0)
end



if(wannaplot)
    figure(1);grid on;
    semilogy(densityBS,pBprime);
    ylim([1e-6,1]);title('Marginal prob of Blockage')
    legend(legendArray);
    
    
    figure(2);grid on;
    semilogy(densityBS,pBCondPrime);
    title('Conditional prob of Bl given n!=0');
    ylim([1e-6,1])
    legend(legendArray);
    figure(3);grid on;
    semilogy(densityBS,freqPrime)
    title('Expected Freq of blockage')
    legend(legendArray);
    ylim([1e-6,1])
    figure(4);grid on;
    semilogy(densityBS,freqCondPrime);
    title('Conditional expectation of freq of bl given n!=0')
    ylim([1e-6,1])
    legend(legendArray);
    
    figure(5); grid on;
    plot(densityBS,durCond2Prime)
    hold on
    %     plot(densityAP',durCond2(:,6),'y-')
    %     plot(densityAP',durCond3(:,6),'b-')
    title('Conditional expectation of duration of bl given n!=0')
    legend(legendArray);
    
end



%% Effect of cell radius
%Just copied the above code with minor modifications

if(wannaplotCellRadius)
    clear
    close all
    wannaplot=1;
    V = 1; %velocity m/s
    hb = 1.8;
    hr = 1.4;
    ht = 5;
    frac = (hb-hr)/(ht-hr); %fraction depends on heights
    mu = 2; %Expected bloc dur =1/mu
    R = 100; %m Radius
    densityBS = [0.01,1,50,100,150,200,250,300,350,400]*10^(-6); %BS=Base Station
    % densityBS = [100]*10^(-6); %BS=Base Station
    densityBL = [ 0.01]; %Dynamic BLockers
    densityD = [0.0001]; %D = static blockage
    omegaVal = [ pi/3];
    % Rvalues = 50:500; %m Radius
    Rvalues = [100, 200];
    nR = length(Rvalues);
    nBS = length(densityBS);
    nBL = length(densityBL);
    nD = length(densityD);
    nO = length(omegaVal);
    El = 10; %m
    Ew = 10; %m
    
    % !!!-------------------------------!!!
    pB = zeros(nBS,nBL,nO);
    pBCond=zeros(nBS,nBL,nO);
    
    freqCond=zeros(nBS,nBL,nO);
    durCond=zeros(nBS,nBL,nO);
    
    for iT = 1:nBS
        tempind = 0;
        for iR = 1:nR
            R = Rvalues(iR);
            
            
            
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
                        a(iB) = 2*mu/(R*C(iB))-2*mu^2/(R^2*C(iB)^2)*log(1+R*C(iB)/mu);
                        
                        % Coverage probability for open park like area
                        pnn0(iT,iO)=exp(-p(iO)*lamT*pi*R^2);
                        
                        beta(iD) = 2/pi*lamD*(El+Ew);
                        beta0(iD) = lamD*El*Ew;
                        
                        b(iD) = 2*exp(-beta0(iD))/(beta(iD)^2*R^2) * ...
                            (1-(1+beta(iD)*R)*exp(-beta(iD)*R));
                        
                        %write code to verify b: verified. Both matches
                        bExp = @(r) exp(-(beta(iD).*r+beta0(iD))) *2.*r/R^2;
                        b(iD) = integral(bExp,0,R);
                        %coverage probability given by static and self blockage
                        pCoverage(iT,iD,iO) = 1-exp(-b(iD)*p(iO)*lamT*pi*R^2);
                        pNoCoverage(iT,iD,iO) = 1-pCoverage(iT,iD,iO);
                        
                        %1. Marginal prob of Blockage (open area)
                        pB(iT,iB,iO) = exp(-a(iB)*p(iO)*lamT*pi*R^2);
                        
                        %1.1 considering static blockage too.
                        c=C(iB)/mu;
                        
                        aPrimeExp = @(r) exp(-(beta(iD)*r+beta0(iD)))./(1+c*r) *2.*r/R^2;
                        aPrime(iB,iD) = integral(aPrimeExp, 0, R);
                        
                        pBprime(iT,tempind) = exp(-aPrime(iB,iD)*p(iO)*lamT*pi*R^2);
                        
                        %2. Conditional prob of Bl given n!=0
                        pBCond(iT,iB,iO) = (pB(iT,iB,iO)-pnn0(iT,iO))./(1-pnn0(iT,iO));
                        
                        %2.1 Conditional prob of blockage given coverage(newly defined)
                        pBCondPrime(iT,tempind) = (pBprime(iT,tempind)-...
                            pNoCoverage(iT,iD,iO))/pCoverage(iT,iD,iO);
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
    
    
    %     writetable(cell2table([colTitle; num2cell([Rvalues', pBCondPrime])]),...
    %         'figures/theory_withR_LOS.csv','writevariablenames',0);
    %       save('figures/R_LOS.mat','pBCondPrime');
    
    
    
    figure(8);grid on;
    semilogy(densityBS,pBCondPrime, 'LineWidth',2);
    xlabel('BS Density');
    title('LOS Blockage Probability');
    ylim([1e-6,1])
    legend(legendArray);
    %     legend('R=100m','R=200m','R=500m','R=1000m');
    
end