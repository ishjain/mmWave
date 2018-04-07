% #plotResults.m
% Take simData.csv and theoryData.csv and plot the results

% allData = [pB,pBgiven,freq,freqCond,durCond];
% It is (5x10)*5 = (5x50); Note the 5 terms above for order. 
% Create 5 plots for those 5 things

% densityBL = [0.01,0.1,0.2,0.5,0.65];
% densityAP = (1:1:10)/10^4;
clear;
close all;
theory_pB=importdata('figures2/theory_pB.csv',',',1);
theory_pBCond=importdata('figures2/theory_pBCond.csv',',',1);
theory_freq=importdata('figures2/theory_freq.csv',',',1);
theory_freqCond=importdata('figures2/theory_freqCond.csv',',',1);
theory_durCond=importdata('figures2/theory_durCond.csv',',',1);

sim_pB=importdata('figures2/sim_pB.csv',',',1);
sim_pBCond=importdata('figures2/sim_pBCond.csv',',',1);
sim_freq=importdata('figures2/sim_freq.csv',',',1);
sim_freqCond=importdata('figures2/sim_freqCond.csv',',',1);
sim_durCond=importdata('figures2/sim_durCond.csv',',',1);

nlines = size(sim_pB.data,2)-1;

color = {'r','g','b','m','k','y'};
% legendArray = {'\lambda_B=0.1, \omega=0','\lambda_B=0.2, \omega=0',...
%     '\lambda_B=0.1, \omega=2\pi/3','\lambda_B=0.2, \omega=2\pi/3'};
legendArray= {'lamB0.01omega0','lamB0.1omega0','lamB0.2omega0',...
    'lamB0.01omega90','lamB0.1omega90','lamB0.2omega90'};
h=figure(1);
grid on;
for j=1:nlines
semilogy(theory_pB.data(:,1), sim_pB.data(:,j+1),'color',color{j},'LineStyle','--','marker','o','LineWidth',2)
hold on;
end
for j=1:nlines
semilogy(theory_pB.data(:,1), theory_pB.data(:,j+1),'color',color{j},'LineStyle','-','marker','*','LineWidth',2)
hold on;
end

ylim([10^-6,1])
g=legend(legendArray);
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
for j=1:nlines
semilogy(theory_pB.data(:,1), sim_pBCond.data(:,j+1),'color',color{j},'LineStyle','--','marker','o','LineWidth',2)
hold on;
end
for j=1:nlines
semilogy(theory_pB.data(:,1), theory_pBCond.data(:,j+1),'color',color{j},'LineStyle','-','marker','*','LineWidth',2)
hold on;
end

ylim([10^-6,1])
g=legend(legendArray);
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
for j=1:nlines
semilogy(theory_pB.data(:,1), sim_freq.data(:,j+1),'color',color{j},'LineStyle','--','marker','o','LineWidth',2)
hold on;
end
for j=1:nlines
semilogy(theory_pB.data(:,1), theory_freq.data(:,j+1),'color',color{j},'LineStyle','-','marker','*','LineWidth',2)
hold on;
end
ylim([10^-6,1])
g=legend(legendArray);
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
for j=1:nlines
semilogy(theory_pB.data(:,1), sim_freqCond.data(:,j+1),'color',color{j},'LineStyle','--','marker','o','LineWidth',2)
hold on;
end
for j=1:nlines
semilogy(theory_pB.data(:,1), theory_freqCond.data(:,j+1),'color',color{j},'LineStyle','-','marker','*','LineWidth',2)
hold on;
end
ylim([10^-6,1])
g=legend(legendArray);
   
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
for j=1:nlines
semilogy(theory_pB.data(:,1), sim_durCond.data(:,j+1),'color',color{j},'LineStyle','--','marker','o','LineWidth',2)
hold on;
end
for j=1:nlines
semilogy(theory_pB.data(:,1), theory_durCond.data(:,j+1),'color',color{j},'LineStyle','-','marker','*','LineWidth',2)
hold on;
end
ylim([10^-2,1])
g=legend(legendArray);
set(g,'fontsize',13)
xlabel('AP Density \lambda_T (x 100/km^2)', 'fontsize',13)
ylabel('Average blockage duration conditioned on n~=0','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures2/durCond.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures2/durCond.png'],'-dpng','-r0')
