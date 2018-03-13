% #plotResults.m
% Take simData.csv and theoryData.csv and plot the results

% allData = [pB,pBgiven,freq,freqCond,durCond];
% It is (5x10)*5 = (5x50); Note the 5 terms above for order. 
% Create 5 plots for those 5 things

densityBL = [0.01,0.1,0.2,0.5,0.65];
densityAP = (1:1:10)/10^4;

theory_pB=csvread('figures2/theory_pB.csv');
theory_pBCond=csvread('figures2/theory_pBCond.csv');
theory_freq=csvread('figures2/theory_freq.csv');
theory_freqCond=csvread('figures2/theory_freqCond.csv');
theory_durCond=csvread('figures2/theory_durCond.csv');

sim_pB=csvread('figures2/sim_pB.csv');
sim_pBCond=csvread('figures2/sim_pBCond.csv');
sim_freq=csvread('figures2/sim_freq.csv');
sim_freqCond=csvread('figures2/sim_freqCond.csv');
sim_durCond=csvread('figures2/sim_durCond.csv');

color = {'r','g','b','m','k'};
h=figure(1);
grid on;
for j=1:5
semilogy(densityAP*10^4, theory_pB(:,j+1),'color',color{j},'LineWidth',2)
hold on;
end
for j=1:5
semilogy(densityAP*10^4, sim_pB(:,j+1),'color',color{j},'LineStyle','--','LineWidth',2)
hold on;
end
ylim([10^-4,1])
g=legend('\rho_B=0.01 bl/m^2','\rho_B=0.1 bl/m^2','\rho_B=0.2 bl/m^2',...
    '\rho_B=0.5 bl/m^2','\rho_B=0.65 bl/m^2' );
set(g,'fontsize',13)
xlabel('AP Density \lambda_T (x 100/km^2)', 'fontsize',13)
ylabel('Prob all APs blocked within 100m','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures2/pB.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures2/pB.png'],'-dpng','-r0')

h=figure(2);
grid on;
for j=1:5
semilogy(densityAP*10^4, theory_pBCond(:,j+1),'color',color{j},'LineWidth',2)
hold on;
end
for j=1:5
semilogy(densityAP*10^4, sim_pBCond(:,j+1),'color',color{j},'LineStyle','--','LineWidth',2)
hold on;
end
ylim([10^-4,1])
g=legend('\rho_B=0.01 bl/m^2','\rho_B=0.1 bl/m^2','\rho_B=0.2 bl/m^2',...
    '\rho_B=0.5 bl/m^2','\rho_B=0.65 bl/m^2' );
set(g,'fontsize',13)
xlabel('AP Density \lambda_T (x 100/km^2)', 'fontsize',13)
ylabel('P(B|n~=0)','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures2/pBCond.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures2/pBCond.png'],'-dpng','-r0')

h=figure(3);
grid on;
for j=1:5
semilogy(densityAP*10^4, sim_freq(:,j+1),'color',color{j},'LineWidth',2)
hold on;
end
for j=1:5
semilogy(densityAP*10^4, theory_freq(:,j+1),'color',color{j},'LineStyle','--','LineWidth',2)
hold on;
end
ylim([10^-4,1])
g=legend('\rho_B=0.01 bl/m^2','\rho_B=0.1 bl/m^2','\rho_B=0.2 bl/m^2',...
    '\rho_B=0.5 bl/m^2','\rho_B=0.65 bl/m^2' );
set(g,'fontsize',13)
xlabel('AP Density \lambda_T (x 100/km^2)', 'fontsize',13)
ylabel('Expected blockage frequency','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures2/freq.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures2/freq'],'-dpng','-r0')

h=figure(4);
grid on;
for j=1:5
semilogy(densityAP*10^4, sim_freqCond(:,j+1),'color',color{j},'LineWidth',2)
hold on;
end
for j=1:5
semilogy(densityAP*10^4, theory_freqCond(:,j+1),'color',color{j},'LineStyle','--','LineWidth',2)
hold on;
end
ylim([10^-4,1])
g=legend('\rho_B=0.01 bl/m^2','\rho_B=0.1 bl/m^2','\rho_B=0.2 bl/m^2',...
    '\rho_B=0.5 bl/m^2','\rho_B=0.65 bl/m^2' );
set(g,'fontsize',13)
xlabel('AP Density \lambda_T (x 100/km^2)', 'fontsize',13)
ylabel('Expected blockage frequency conditioned on n~=0','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures2/freq.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures2/freq'],'-dpng','-r0')

h=figure(5);
grid on;
for j=1:5
semilogy(densityAP*10^4, sim_durCond(:,j+1),'color',color{j},'LineWidth',2)
hold on;
end
for j=1:5
semilogy(densityAP*10^4, theory_durCond(:,j+1),'color',color{j},'LineStyle','--','LineWidth',2)
hold on;
end
ylim([10^-3,1])
g=legend('\rho_B=0.01 bl/m^2','\rho_B=0.1 bl/m^2','\rho_B=0.2 bl/m^2',...
    '\rho_B=0.5 bl/m^2','\rho_B=0.65 bl/m^2' );
set(g,'fontsize',13)
xlabel('AP Density \lambda_T (x 100/km^2)', 'fontsize',13)
ylabel('Average blockage duration conditioned on n~=0','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures2/durCond.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures2/durCond.png'],'-dpng','-r0')
