%casestudy.m
% Mar 19: we perform case study of AR/VR application
% Plot AP density and AP height with Bl density for 3 case
% 1. Bl probability < p_bar
% 2. Bl freq < freq_bar
% 3. Bl duration < dur_bar

clear
close all

wannaplot=1;
V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 5;
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
densityBL = 0.01:0.01:.20;
densityAP = 10:0.1:500*1e-6;%[10,20,50,100,200,500]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, 2*pi/3];

pB_bar = 1e-5;%0.001; % probbility
freq_bar = 1/(60*60*24);%1/60;%1/200; %intrrutions/sec
% freqCond_bar = 0.001;
dur_bar = 0.040; %sec
omega=0;
for indB = 1:length(densityBL)
C = 2/pi*densityBL(indB)*V*frac;


a_old = 1-2*mu./(R*C)+2*mu^2./(R^2*C.^2).*log(1+R*C/mu);
p = omega/(2*pi);
a= a_old*(1-p)+p;
% lamT = log(p_bar)./(a-1)/(pi*R^2);


%Note I removed R^2 so change lamT*pi*R^2 to lamT*pi (only change in unit)
pB = @(lamT) exp((a-1).*lamT.*pi)-pB_bar;
pBCond = @(lamT) (exp((a-1).*lamT.*pi)-exp(-lamT.*pi))./(1-exp(-lamT.*pi))-pB_bar;
freq = @(lamT) mu*a*lamT*pi*exp((a-1)*lamT*pi) - freq_bar;
freqCond = @(lamT) (mu*a*lamT*pi*exp((a-1)*lamT*pi))./(1-exp(-lamT.*pi)) - freq_bar;
durCond = @(lamT) exp(-lamT.*pi)*(ei(lamT*pi)-log(lamT*pi)-0.5772)/(mu*(1-exp(-lamT.*pi)))-dur_bar;

% pB = @(lamT) fun_pB(lamT);
% pBCond = @(lamT) fun_pBCond(lamT,a,pB_bar);
lamTxx(indB) = -log(pB_bar).*(1+2*R*C/(3*mu))/(pi*R^2); %approx
lamT(indB) = fzero(pB,[0.01,100]);
lamT2(indB) = fzero(pBCond,[0.01,20]);
% lamT_freq(indB) = fzero(freq,[0.00001,200000]);
lamT_freqCond(indB) = fzero(freqCond,[0.01,20]);
lamT_durCond(indB) = fzero(durCond,[0.1,200]);
end


% figure(1);hold on;grid on;
%  plot(densityBL,lamT,'linewidth',1);
%  plot(densityBL, lamTxx,'--')
%  xlabel('Blocker Density (\lambda_T bl/m^2)')
% ylabel('BS Density (\lambda_T) x100/km^2')
%  legend('Theory (P(B)=1e-5)', 'Linear approx (P(B)=1e-5)')
h= figure(3); hold on; grid on
plot(densityBL,lamT2,'linewidth',2);
xlabel('Blocker density (\lambda_B bl/m^2)','fontsize',13)
ylabel('BS Density (\lambda_T) x100/km^2','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/BS_vs_bl_density.pdf'],'-dpdf','-r0')
print(h,[pwd '/BS_vs_bl_density.png'],'-dpng','-r0')