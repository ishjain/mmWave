%BS_height_density_tradeoff.m
%Apr 8: For lamB = [0.01, 0.1, 0.2], plot the h_T vs lamT tradeoff

clear
close all

wannaplot=1;
V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
htVal = 2:10;

mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
densityBL = [0.1,0.01];
% densityAP = 10:0.1:500*1e-6;%[10,20,50,100,200,500]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, pi/3];

pB_bar = 1e-5;%0.001; % probbility
freq_bar = 1e-7;%1/60;%1/200; %intrrutions/sec
% freqCond_bar = 0.001;
dur_bar = 0.020; %sec
for indO=1:length(omegaVal)
    for indB = 1:length(densityBL)
        for indH = 1:length(htVal)
            omega = omegaVal(indO);
            frac(indH) = (hb-hr)/(htVal(indH)-hr);
            C = 2/pi*densityBL(indB)*V*frac(indH);
            
            
            a = 2*mu./(R*C)-2*mu^2./(R^2*C.^2).*log(1+R*C/mu);
            p = 1-omega/(2*pi);
            % a= a_old*(1-p)+p;
            % lamT = log(p_bar)./(a-1)/(pi*R^2);
            
            
            %Note I removed R^2 so change lamT*pi*R^2 to lamT*pi (only change in unit)
            pB = @(lamT) exp(-a*p*lamT*pi)-pB_bar;
            pBCond = @(lamT) (exp(-a*p*lamT.*pi)-exp(-p*lamT.*pi))./(1-exp(-p*lamT*pi))-pB_bar;
            freq = @(lamT) mu*(1-a)*p*lamT*pi*exp(-a*p*lamT*pi) - freq_bar;
            freqCond = @(lamT) (mu*(1-a)*p*lamT*pi*exp(-a*p*lamT*pi))./(1-exp(-p*lamT*pi)) - freq_bar;
            durCond = @(lamT) exp(-p*lamT*pi)*(ei(p*lamT*pi)-log(p*lamT*pi)-0.5772)/(mu*(1-exp(-p*lamT.*pi)))-dur_bar;
            
            % pB = @(lamT) fun_pB(lamT);
            % pBCond = @(lamT) fun_pBCond(lamT,a,pB_bar);
            lamTxx(indO,indB,indH) = -log(pB_bar).*(1+2*R*C/(3*mu))/(pi*R^2); %approx
            lamT_pB(indO,indB,indH) = fzero(pB,[0.01,100]);
            lamT_pBCond(indO,indB,indH) = fzero(pBCond,[0.01,200]);
            % lamT_freq(indO,indB,indH) = fzero(freq,[0.00001,200000]);
            lamT_freqCond(indO,indB,indH) = fzero(freqCond,[0.01,200]);
            lamT_durCond(indO,indB,indH) = fzero(durCond,[0.1,200]);
        end
    end
end
lamT_pBCond=reshape(lamT_pBCond,4,length(htVal));
% figure(1);hold on;grid on;
% plot(densityBL,lamT_pB,'linewidth',1);
% plot(densityBL, lamTxx,'--')
% xlabel('Blocker Density (\lambda_T bl/m^2)')
% ylabel('BS Density (\lambda_T) x100/km^2')
% legend('Theory (P(B)=1e-5)', 'Linear approx (P(B)=1e-5)')

xVar = htVal;



colTitle= {'BSheight','lambdaB0.1omega0','lambdaB0.1omega60',...
    'lambdaB0.01omega0','lambdaB0.01omega60'};

%%Save to csv file for Rajeev
writetable(cell2table([colTitle; num2cell([xVar',lamT_pBCond'])]),...
    [pwd '\heightvsdensity.csv'],'writevariablenames',0)


colorVal = {'m','c','r','b'};
markerVal = {'o','d','^','s'};
h=figure(2); hold on; grid on;
% plot(xVar,lamT_pBCond(1,:),'rd-','linewidth',2)
% plot(xVar,lamT_pBCond(3,:),'bs-','linewidth',2)
for ii=1:length(densityBL)*length(omegaVal)
    plot(xVar,lamT_pBCond(ii,:),'linewidth',2,'color',colorVal{ii},...
        'marker',markerVal{ii});
end
g=legend('\lambda_B=0.1, \omega=0','\lambda_B=0.1, \omega=\pi/3',...
    '\lambda_B=0.01, \omega=0','\lambda_B=0.01, \omega=\pi/3');
% g=legend('\lambda_B=0.1','\lambda_B=0.01');
set(g,'fontsize',13)
xlabel('BS height (h_T) m','fontsize',13)
ylabel('BS Density (\lambda_T) x100/km^2','fontsize',13)
% title('BS height-density tradeoff for P(B|C) = 1e-7')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/heightvsdensity.pdf'],'-dpdf','-r0')
print(h,[pwd '/heightvsdensity.png'],'-dpng','-r0')