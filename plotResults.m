% #plotResults.m
% Take simData.csv and theoryData.csv and plot the results

% allData = [pB,pBgiven,freq,freqCond,durCond];
% It is (5x10)*5 = (5x50); Note the 5 terms above for order. 
% Create 5 plots for those 5 things

densityBL = [0.01,0.1,0.2,0.5,0.65];
densityAP = (1:1:10)/10^4;
simData = csvread('simData.csv');


color = {'r','g','b','m','k'};
h=figure(1);
grid on;
for j=1:5
semilogy(densityAP*10^4, probAllBl(j,:),'color',color{j},'LineWidth',2)
hold on;
end
for j=1:5
semilogy(densityAP*10^4, th_probAllBl(j,:),'color',color{j},'LineStyle','--','LineWidth',2)
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
print(h,[pwd '/figures2/prob_all_bl_new.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures2/prob_all_bl_new.png'],'-dpng','-r0')

h=figure(2);
grid on;
for j=1:5
semilogy(densityAP*10^4, avgFreq(j,:),'color',color{j},'LineWidth',2)
hold on;
end
for j=1:5
semilogy(densityAP*10^4, th_freqBl(j,:),'color',color{j},'LineStyle','--','LineWidth',2)
hold on;
end
ylim([10^-4,1])
g=legend('\rho_B=0.01 bl/m^2','\rho_B=0.1 bl/m^2','\rho_B=0.2 bl/m^2',...
    '\rho_B=0.5 bl/m^2','\rho_B=0.65 bl/m^2' );
set(g,'fontsize',13)
xlabel('AP Density \lambda_T (x 100/km^2)', 'fontsize',13)
ylabel('Average frequency of all APs blockage','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures2/avg_freq_new.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures2/avg_freq_new.png'],'-dpng','-r0')

h=figure(3);
grid on;
for j=1:5
semilogy(densityAP*10^4, avgDur(j,:),'color',color{j},'LineWidth',2)
hold on;
end
for j=1:5
semilogy(densityAP*10^4, th_durBl(j,:),'color',color{j},'LineStyle','--','LineWidth',2)
hold on;
end
ylim([10^-3,1])
g=legend('\rho_B=0.01 bl/m^2','\rho_B=0.1 bl/m^2','\rho_B=0.2 bl/m^2',...
    '\rho_B=0.5 bl/m^2','\rho_B=0.65 bl/m^2' );
set(g,'fontsize',13)
xlabel('AP Density \lambda_T (x 100/km^2)', 'fontsize',13)
ylabel('Average blockage duration of all APs','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures2/avg_dur_new.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures2/avg_dur_new.png'],'-dpng','-r0')
