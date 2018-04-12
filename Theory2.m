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
frac = (hb-hr)/(ht-hr);
% simTime = 600*10; %sec Total Simulation time
% tstep = 0.0001; %(sec) time step
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
temp = 2/pi*V*frac*R;
0^4;
densityBL = [0.01,0.1,0.2];
densityAP = [50,100,200,300,400,500,600]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, pi/3, pi/2];


%calculat Ei by wolform alpha for given lambdaAP and R=100;
% Ei = [9.21,105.059, 1499.93, 25031.9, 453682, 8.6341e6, 1.697e8, 3.41357e9, 6.9864e10, 1.4493e12];

% !!!-------------------------------!!!
% pB = zeros(length(densityAP),length(densityBL));
% pBCond=zeros(length(densityAP),length(densityBL));
% freq=zeros(length(densityAP),length(densityBL));
% freqCond=zeros(length(densityAP),length(densityBL));
% durCond=zeros(length(densityAP),length(densityBL));

for indT = 1:length(densityAP)
    for indB = 1:length(densityBL)
        for indO = 1:length(omegaVal)
            lamT = densityAP(indT); %lambda
            lamB = densityBL(indB);
            omega = omegaVal(indO);
            p = omega/(2*pi);
%             angFrac = 1-omega/(2*pi);
            C(indB) = 2/pi*lamB*V*frac;
            a(indB) = 1-2*mu/(R*C(indB))+2*mu^2/(R^2*C(indB)^2)*log(1+R*C(indB)/mu);
            pnn0(indT,indO)=exp(-(1-p)*lamT.*pi*R^2);%*angFrac);
            
            %for self blockage, 
            
            a_tilde(indB,indO) = a(indB)*(1-p)+p;
            %1. Marginal prob of Blockage
            pB(indT,indB,indO) = exp((a_tilde(indB,indO)-1).*lamT.*pi*R^2);%*angFrac;% + 1-angFrac;
            pB_tilde(indT,indB,indO) = exp((a_tilde(indB,indO)-1).*lamT.*pi*R^2);
            %2. Conditional prob of Bl given n!=0
            pBCond(indT,indB,indO) = (pB(indT,indB,indO)-pnn0(indT,indO))./(1-pnn0(indT,indO));
            
            %3. Expected Freq of blockage
            freq(indT,indB,indO) = mu*a_tilde(indB,indO)*lamT*pi*R^2*exp((a_tilde(indB,indO)-1)*lamT*pi*R^2);
            
            %4. Conditional expectation of freq of bl given n!=0
            freqCond(indT,indB,indO) = (freq(indT,indB,indO))./(1-pnn0(indT,indO));
            
            %5. Conditional expectation of duration of bl given n!=0
            
            
            %This is used to get
            b(indT,indO) = lamT*pi*R^2;%*angFrac;
%             ab(indT,indB,indO) = b(indT,indO)*a_tilde(indB,indO);%*b(indT);
            Ei(indT,indB,indO) = ei(b(indT,indO))-log(b(indT,indO))-0.5772;
            
            durCond(indT,indB,indO) = pnn0(indT,indO)*Ei(indT,indB,indO)/(mu*(1-pnn0(indT,indO)));
            %approx
            durCond2(indT,indB,indO) = (1/b(indT,indO)+3/b(indT,indO)^2)/...
                (mu*(1-pnn0(indT,indO)));
            durCond2(indT,indB,indO) = (1/b(indT,indO))/...
                (mu*(1-pnn0(indT,indO)));
            
            %         nn = poissrnd(lamT*pi*R^2*angFrac,1,1000);
            %         nn_n0 = nn(nn~=0);
            %         durCond2(indT,indB,indO) = mean(1./(nn_n0*2));%verified
            %         durCond3(indT,indB,indO) = pBCond(indT,indB,indO)/freqCond(indT,indB,indO);%wrong
        end
    end
end




%save csv file
%%Make it 2D arrays
pB=reshape(pB(:,:,[1,3]),length(densityAP),6);
pB_tilde=reshape(pB_tilde(:,:,[1,3]),length(densityAP),6);
pBCond=reshape(pBCond(:,:,[1,3]),length(densityAP),6);
freq=reshape(freq(:,:,[1,3]),length(densityAP),6);
freqCond=reshape(freqCond(:,:,[1,3]),length(densityAP),6);
durCond=reshape(durCond(:,:,[1,3]),length(densityAP),6);

%6 rows represents 6 values of AP density
%8 columns has the following legend
% legendArray = {'\lambda_B=0.1,\omega=0','\lambda_B=0.2,\omega=0',...
%     '\lambda_B=0.1,\omega=\pi/3','\lambda_B=0.2,\omega=\pi/3',...
%     '\lambda_B=0.1,\omega=\pi/2','\lambda_B=0.2,\omega=\pi/2',...
%     '\lambda_B=0.1,\omega=2\pi/3','\lambda_B=0.2,\omega=2\pi/3'};
% legendArray = {'\lambda_B=0.1, \omega=0','\lambda_B=0.2, \omega=0',...
%     '\lambda_B=0.1, \omega=2\pi/3','\lambda_B=0.2, \omega=2\pi/3'};
legendArray= {'lamB0.01omega0','lamB0.1omega0','lamB0.2omega0',...
    'lamB0.01omega90','lamB0.1omega90','lamB0.2omega90'};
colTitle= {'lamT','lamB0.01omega0','lamB0.1omega0','lamB0.2omega0',...
    'lamB0.01omega90','lamB0.1omega90','lamB0.2omega90'};
%     'lamB0.1omega2pi3','lamB0.2omega2pi3'};

%%Save to csv file for Rajeev
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,pB])]),...
    'figures2/theory_pB.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,pBCond])]),...
    'figures2/theory_pBCond.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,freq])]),...
    'figures2/theory_freq.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,freqCond])]),...
    'figures2/theory_freqCond.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,durCond])]),...
    'figures2/theory_durCond.csv','writevariablenames',0)

% csvwrite('figures2/theory_pB.csv',[densityAP'*10^4,pB]);
% csvwrite('figures2/theory_pBCond.csv',[densityAP'*10^4,pBCond]);
% csvwrite('figures2/theory_freq.csv',[densityAP'*10^4,freq]);
% csvwrite('figures2/theory_freqCond.csv',[densityAP'*10^4,freqCond]);
% csvwrite('figures2/theory_durCond.csv',[densityAP'*10^4,durCond]);


if(wannaplot)
    figure(1);
    semilogy(densityAP,pB);
    ylim([1e-4,1]);title('Marginal prob of Blockage')
    legend(legendArray);
    figure(2);
    semilogy(densityAP,pBCond); title('Conditional prob of Bl given n!=0')
    ylim([1e-4,1])
    legend(legendArray);
    figure(3);
    semilogy(densityAP,freq)
    title('Expected Freq of blockage')
    legend(legendArray);
    ylim([1e-4,1])
    figure(4);
    semilogy(densityAP,freqCond);
    title('Conditional expectation of freq of bl given n!=0')
    ylim([1e-4,1])
    legend(legendArray);
    figure(6); grid on;
    plot(densityAP,durCond(:,1:2))
    hold on
    plot(densityAP(3:6),durCond2(3:6,1:2))
    title('Conditional expectation of duration of bl given n!=0')
    legend(legendArray);
    
end




