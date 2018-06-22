%TheoryLOS: June19: Use generalized model and find dynamic bl probab given
%coverage which is defined as "atleast one BS out side blocked zone
%(defined by static ans self-blockage)"
%Theory_LOS: June6: Renamed Theory2 to Theory_LOS and tehn to TheoryLOS
%Theory2 : Mar16: Based on new range of densities and self blockage angle
%Based on Data7 onweards
%Mar10: Theoretically calculating various values and save csv file

clear
close all
%Copy from Simulation.m !!!---------!!!
wannaplot=1;
V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 5;
frac = (hb-hr)/(ht-hr); %fraction depends on heights
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
densityBS = [50,100,200,300,400]*10^(-6); %BS=Base Station
densityBL = [0.01, 0.1]; %Dynamic BLockers
densityD = [1e-9,0.0001]; %D = static blockage
omegaVal = [0 pi/3];

nBS = length(densityBS);
nBL = length(densityBL);
nD = length(densityD);
nO = length(omegaVal);
El = 10; %m
Ew = 10; %m

% !!!-------------------------------!!!
pB = zeros(nBS,nBL,nO);
pBCond=zeros(nBS,nBL,nO);
freq=zeros(nBS,nBL,nO);
freqCond=zeros(nBS,nBL,nO);
durCond=zeros(nBS,nBL,nO);

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
                %                 aPrime(iB,iD) = 2*mu*exp(-beta0(iD))/R/C(iB)*...
                %                     (mu*exp(beta(iD)*mu/C(iB))/R/C(iB)*(ei(-beta(iD)/c)-ei(beta(iD)*R-beta(iD)/c))+...
%                     (1-exp(-beta(iD)*R))/R/beta(iD));
                pBprime(iT,tempind) = exp(-aPrime(iB,iD)*p(iO)*lamT*pi*R^2);
                
                %2. Conditional prob of Bl given n!=0
                pBCond(iT,iB,iO) = (pB(iT,iB,iO)-pnn0(iT,iO))./(1-pnn0(iT,iO));
                
                %2.1 Conditional prob of blockage given coverage(newly defined)
                pBCondPrime(iT,tempind) = (pBprime(iT,tempind)-...
                    pNoCoverage(iT,iD,iO))/pCoverage(iT,iD,iO);
                
                %3. Expected Freq of blockage
                freq(iT,iB,iO) = mu*(1-a(iB))*p(iO)*lamT*pi*R^2*exp(-a(iB)*p(iO)*lamT*pi*R^2);
                freqPrime(iT,tempind) = mu*(1-aPrime(iB,iD))*p(iO)*b(iD)...
                    *lamT*pi*R^2*exp(-aPrime(iB,iD)*p(iO)*b(iD)*lamT*pi*R^2);
                %4. Conditional expectation of freq of bl given n!=0
                freqCond(iT,iB,iO) = (freq(iT,iB,iO))./(1-pnn0(iT,iO));
                freqCondPrime(iT,tempind) = freqPrime(iT,tempind)/pCoverage(iT,iD,iO);
                %5. Conditional expectation of duration of bl given n!=0
                
                
                %This is used to get
                tempb(iT,iD,iO) = p(iO)*b(iD)*lamT*pi*R^2;%*angFrac;
                %             ab(indT,indB,indO) = b(indT,indO)*a_tilde(indB,indO);%*b(indT);
                Ei(iT,iB,iO) = ei(tempb(iT,iD,iO))-log(tempb(iT,iD,iO))-0.5772;
%                  durCond(indT,indB,indO) = pnn0(indT,indO)*Ei(indT,indB,indO)/(mu*(1-pnn0(indT,indO)));
                durCondPrime(iT,tempind) = (1-pCoverage(iT,iD,iO))*Ei(iT,iB,iO)/(mu*pCoverage(iT,iD,iO));
                durCondPrime(iT,tempind)=durCondPrime(iT,tempind)*1000;
                %approx
                durCond(iT,iB,iO) = (1/tempb(iT,iD,iO)+1/tempb(iT,iD,iO)^2)/...
                    (mu*(1-pnn0(iT,iO)));
                durCond2(iT,iB,iO) = (1/tempb(iT,iD,iO))/...
                    (mu*(1-pnn0(iT,iO)));
                durCond2prime(iT,iB,iO) = (1/tempb(iT,iD,iO))/...
                    (mu*(1-pnn0(iT,iO)));
                
                %         nn = poissrnd(lamT*pi*R^2*angFrac,1,1000);
                %         nn_n0 = nn(nn~=0);
                %         durCond2(indT,indB,indO) = mean(1./(nn_n0*2));%verified
                %         durCond3(indT,indB,indO) = pBCond(indT,indB,indO)/freqCond(indT,indB,indO);%wrong
                colTitle{1}='lamT';
                if(lamD<1e-7),lamD=0;end
%                 colTitle{tempind+1} = sprintf('St%dDy%dOm%d',iB,iD,iO);
                colTitle{tempind+1} = strcat('lamB',num2str(lamB),...
                    'lamD',num2str(lamD*1e4),'omega',num2str(omega*360/2/pi));
                legendArray{tempind} = strcat(' \lambda_B=',num2str(lamB),...
                    '\lambda_D=',num2str(lamD*1e4),' \omega=',num2str(omega*360/2/pi));
%                 sprintf('\lambda_B=%f,\lambda_D=%d,\omega=%d',lamB,lamD*1e6,omega*360/2/pi);
            end
        end
    end
end



%save csv file
%%Make it 2D arrays
% pB=reshape(pB(:,:,:),nBS,nBL*nO);
% % pB_tilde=reshape(pB_tilde(:,:,:),nBS,nBL*nO);
% pBCond=reshape(pBCond(:,:,:),nBS,nBL*nO);
% freq=reshape(freq(:,:,:),nBS,nBL*nO);
% freqCond=reshape(freqCond(:,:,:),nBS,nBL*nO);
% durCond=reshape(durCond(:,:,:),nBS,nBL*nO);
% durCond2=reshape(durCond2(:,:,:),nBS,nBL*nO);
% durCond3=reshape(durCond3(:,:,:),nBS,nBL*nO);

% pBprime_ = reshape(pBprime,nBS,nBL*nD*nO);
% pBCondPrime_ = reshape(pBCondPrime,nBS,nBL*nD*nO);
% freqPrime_ = reshape(freqPrime,nBS,nBL*nD*nO);
% freqCondPrime_ = reshape(freqCondPrime,nBS,nBL*nD*nO);
% durCondPrime_ = reshape(durCondPrime,nBS,nBL*nD*nO);
% legendArray = {'\lambda_B=0.1,\omega=0','\lambda_B=0.2,\omega=0',...
%     '\lambda_B=0.1,\omega=\pi/3','\lambda_B=0.2,\omega=\pi/3',...
%     '\lambda_B=0.1,\omega=\pi/2','\lambda_B=0.2,\omega=\pi/2',...
%     '\lambda_B=0.1,\omega=2\pi/3','\lambda_B=0.2,\omega=2\pi/3'};
% legendArray = {'\lambda_B=0.1, \omega=0','\lambda_B=0.2, \omega=0',...
%     '\lambda_B=0.1, \omega=2\pi/3','\lambda_B=0.2, \omega=2\pi/3'};
% legendArray= {'lamB0.01omega0','lamB0.1omega0','lamB0.2omega0',...
%     'lamB0.01omega90','lamB0.1omega90','lamB0.2omega90'};
% colTitle= {'lamT','lamB0.01omega0','lamB0.1omega0','lamB0.2omega0',...
%     'lamB0.01omega60','lamB0.1omega60','lamB0.2omega60'};
%     'lamB0.1omega2pi3','lamB0.2omega2pi3'};

% colTitle= {'lamT','lamB0.01omega0','lamB0.1omega0',...
%     'lamB0.01omega60','lamB0.1omega60'};
% legendArray = {'lamB0.01omega0','lamB0.1omega0',...
%     'lamB0.01omega60','lamB0.1omega60'};
% legendArray = {'lamB0.01omega0','lamB0.1omega0',...
%     'lamB0.01omega60','lamB0.1omega60'};
%Save to csv file for Rajeev
writetable(cell2table([colTitle; num2cell([densityBS'*10^4,pBprime])]),...
    'figures2/theory_pB.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityBS'*10^4,pBCondPrime])]),...
    'figures2/theory_pBCond.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityBS'*10^4,freqPrime])]),...
    'figures2/theory_freq.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityBS'*10^4,freqCondPrime])]),...
    'figures2/theory_freqCond.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityBS'*10^4,durCondPrime])]),...
    'figures2/theory_durCond.csv','writevariablenames',0)

% csvwrite('figures2/theory_pB.csv',[densityAP'*10^4,pB]);
% csvwrite('figures2/theory_pBCond.csv',[densityAP'*10^4,pBCond]);
% csvwrite('figures2/theory_freq.csv',[densityAP'*10^4,freq]);
% csvwrite('figures2/theory_freqCond.csv',[densityAP'*10^4,freqCond]);
% csvwrite('figures2/theory_durCond.csv',[densityAP'*10^4,durCond]);


if(wannaplot)
    figure(1);
    semilogy(densityBS,pBprime);
    ylim([1e-4,1]);title('Marginal prob of Blockage')
    legend(legendArray);
    
    
    figure(2);
    semilogy(densityBS,pBCondPrime); 
    title('Conditional prob of Bl given n!=0');
    ylim([1e-4,1])
    legend(legendArray);
    figure(3);
    semilogy(densityBS,freqPrime)
    title('Expected Freq of blockage')
    legend(legendArray);
    ylim([1e-4,1])
    figure(4);
    semilogy(densityBS,freqCondPrime);
    title('Conditional expectation of freq of bl given n!=0')
    ylim([1e-4,1])
    legend(legendArray);
    
    figure(5); grid on;
    plot(densityBS,durCondPrime)
    hold on
    %     plot(densityAP',durCond2(:,6),'y-')
    %     plot(densityAP',durCond3(:,6),'b-')
    title('Conditional expectation of duration of bl given n!=0')
    legend(legendArray);
    
%      figure(6); grid on;
%     semilogy(densityBS,pBprime_)
%     hold on
%     title('Conditional probability')
% %     legend(legendArray);
% 
%     
%     figure(7); grid on;
%     semilogy(densityBS,pBCondPrime_)
%     hold on
%     title('probability blocage')
end




