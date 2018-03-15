%Theory
%Mar10: Theoretically calculating various values and save csv file

clear
close all
%Copy from Simulation.m !!!---------!!!
wannaplot=1;
V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 6;
frac = (hb-hr)/(ht-hr);
simTime = 60*10; %sec Total Simulation time
tstep = 0.0001; %(sec) time step
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
temp = 2/pi*V*frac*R;
densityBL = [0.01,0.1,0.2,0.5,0.65];
densityAP = (1:1:10)/10^4;
 
%calculat Ei by wolform alpha for given lambdaAP and R=100;
% Ei = [9.21,105.059, 1499.93, 25031.9, 453682, 8.6341e6, 1.697e8, 3.41357e9, 6.9864e10, 1.4493e12];

% !!!-------------------------------!!!
pB = zeros(length(densityAP),length(densityBL));
pBgiven=zeros(length(densityAP),length(densityBL));
freq=zeros(length(densityAP),length(densityBL));
freqCond=zeros(length(densityAP),length(densityBL));
durCond=zeros(length(densityAP),length(densityBL));

for indT = 1:length(densityAP)
    for indB = 1:length(densityBL)
        lamT = densityAP(indT); %lambda
        lamB = densityBL(indB); 
        C(indB) = 2/pi*lamB*V*frac;
        a(indB) = 1-2*mu/(R*C(indB))+2*mu^2/(R^2*C(indB)^2)*log(1+R*C(indB)/mu);
        pnn0(indT)=exp(-lamT.*pi*R^2);
        %1. Marginal prob of Blockage
        pB(indT,indB) = exp((a(indB)-1).*lamT.*pi*R^2);
        
        %2. Conditional prob of Bl given n!=0
        pBgiven(indT,indB) = (pB(indT,indB)-pnn0(indT))./(1-pnn0(indT));
        
        %3. Expected Freq of blockage
        freq(indT,indB) = mu*a(indB)*lamT*pi*R^2*exp((a(indB)-1)*lamT*pi*R^2);
        
        %4. Conditional expectation of freq of bl given n!=0
        freqCond(indT,indB) = (freq(indT,indB))./(1-pnn0(indT));
        
        %5. Conditional expectation of duration of bl given n!=0
        

        %This is used to get 
        b(indT) = lamT*pi*R^2;
        ab(indT,indB) = b(indT);% a(indB)*b(indT);
        Ei(indT,indB) = ei(ab(indT,indB))-log(ab(indT,indB))-0.5772;
        
        durCond1(indT,indB) = pnn0(indT)*Ei(indT,indB)/(mu*(1-pnn0(indT)));
        
        nn = poissrnd(lamT*pi*R^2,1,1000);
        nn_n0 = nn(nn~=0);
        durCond(indT,indB) = mean(1./(nn_n0*2));
        durCond3(indT,indB) = pBgiven(indT,indB)/freqCond(indT,indB);
    end
end

%save csv file
allData = [pB,pBgiven,freq,freqCond,durCond];
csvwrite('figures2/theory_pB.csv',[densityAP;pB']');
csvwrite('figures2/theory_pBCond.csv',[densityAP;pBgiven']');
csvwrite('figures2/theory_freq.csv',[densityAP;freq']');
csvwrite('figures2/theory_freqCond.csv',[densityAP;freqCond']');
csvwrite('figures2/theory_durCond.csv',[densityAP;durCond']');
%Plot
if(wannaplot)
    figure(1);
    semilogy(densityAP,pB); 
    ylim([1e-4,1]);title('Marginal prob of Blockage')
    
    figure(2);
    semilogy(densityAP,pBgiven); title('Conditional prob of Bl given n!=0')
    ylim([1e-4,1])
    
    figure(3);
    semilogy(densityAP,freq)
    title('Expected Freq of blockage')
%     ylim([1e-4,1])
        figure(4);
    semilogy(densityAP,freqCond);
    title('Conditional expectation of freq of bl given n!=0')
%     ylim([1e-4,1])
        figure(5);
   semilogy(densityAP,durCond)
   hold on
    semilogy(densityAP,durCond1)
    semilogy(densityAP,durCond3)
    legend('new','old','divide')
    title('Conditional expectation of duration of bl given n!=0')
    
    
end







