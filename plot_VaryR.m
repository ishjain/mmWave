%plot_VaryR

clear;
close all;
%%Copy from TheoryLOS and THeoryNLOS2

Rvalues = [100, 200, 500, 1000];
densityBS = [10,50,100,200,300,400]*10^(-2); %density/10000
dataNLOS = load('figures2/R_NLOS.mat');
dataLOS = load('figures2/R_LOS.mat');


colorVal = {'r','g','b','k'};
h=figure(1);
grid on;
for ii=1:length(Rvalues)
    semilogy(densityBS,dataLOS.pBCondPrime(ii,:),'color',colorVal{ii}, 'LineWidth',2); hold on;
end

% h= semilogy(densityBS,dataLOS.pBCondPrime, 'LineWidth',2); 
% hold on;
% semilogy(densityBS,dataNLOS.pBCondMin,'--','LineWidth',2);
for ii=1:length(Rvalues)
   semilogy(densityBS,dataNLOS.pBCondMin(ii,:),'--', 'color',colorVal{ii},'LineWidth',2); hold on;
end
    xlabel('BS Density (x100/km^2)','fontsize',14);
    ylabel('Blockage Probability','fontsize',14);
    ylim([1e-6,1])
    legend(['R=',num2str(Rvalues(1)),'m'],['R=',num2str(Rvalues(2)),'m'],...
        ['R=',num2str(Rvalues(3)),'m'],['R=',num2str(Rvalues(4)),'m']); 
set(gca,'fontsize',12)
grid on;

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures2/varyR.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures2/varyR.png'],'-dpng','-r0')

