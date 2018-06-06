% temp
clear
close all
% a vs RC/mu

V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 5;
frac = (hb-hr)/(ht-hr);
% simTime = 600*10; %sec Total Simulation time
% tstep = 0.0001; %(sec) time step
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
omega=0;
densityBL = 0.01:0.01:0.1;
C = 2/pi*densityBL*V*frac;
rate = C([1,10])*100;
a = 2*mu./(R*C)-2*mu^2./(R^2*C.^2).*log(1+R*C/mu);
p = 1-omega/(2*pi);

h=figure(2); hold on; grid on;
plot(densityBL,a,'rd-','linewidth',2)

% g=legend('\lambda_B=0.1','\lambda_B=0.01');
% set(g,'fontsize',13)
xlabel('Blocker densoty (\lambda_B) bl/m^2','fontsize',13)
ylabel('a','fontsize',13)
% title('BS height-density tradeoff for P(B|C) = 1e-7')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/a_vs_density.pdf'],'-dpdf','-r0')
print(h,[pwd '/a_vs_density.png'],'-dpng','-r0')