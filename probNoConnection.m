%probNoConnection
close all
clear
R=100;
rhoT = (.1:.1:2)/100^2;
prob = exp(-rhoT*pi*R^2);
% R=sqrt(10^6/pi);
h=figure(2);
% plot(Lval,emp_noapprox,'LineWidth',2);
hold on; grid on;
plot(rhoT*10^6,prob,'LineWidth',2)
% legend('emperical', 'theoretical')
xlabel('\lambda_T (Density of APs per km^2)', 'fontsize',13)
ylabel('Probability of no connection within 100 m','fontsize',13)
% title('APs blockage probability in 1km^2 area')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures/probNoConnection.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures/probNoConnection.png'],'-dpng','-r0')
% saveas(gcf,[pwd '/figures/probNoConnection.png'])