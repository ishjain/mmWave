close all;
clear all;
v = 1; %m/s
dm = 0.5; %m
EL = dm;
ET = EL/v;
hB = 1.7;
hT = 5;
hR = 1.3;
wS=50; %m
figure;hold on; grid on;
r0val=[50,100];
for ii=1:2
    r0=r0val(ii);
r = r0*(hB-hR)/(hT-hR)+dm/2;
wE = r; %because alpha=0
lambdaI = 1:.05:2.5;
lambda = wE/wS*lambdaI;
Eeta = 1./lambda.*(exp(lambda*ET)-1);
alpha = 1.9; %exponential blockage rate
Eeta2 = 1/alpha;

Eomega = 1./lambda;
Eblockage = 100*Eeta./(Eeta+Eomega);
plot(lambdaI,Eblockage,'LineWidth',2)

Eblockage2 = 100*Eeta2./(Eeta2+Eomega);
plot(lambdaI,Eblockage2,'--','LineWidth',2)
end

leg=legend('r_0=50 M/GI/\infty process',strcat('r_0=50 exp blockage (rate=',num2str(alpha),')'),...
    'r_0=100 M/GI/\infty process', strcat('r_0=100 exp blockage (rate=',num2str(alpha),')')');
set(leg,'FontSize',11,'location','northwest')
xlabel('\lambda_I (Blocker''s arrival rate)','FontSize',13)
ylabel('Fraction of Blockage (%)','FontSize', 13)
saveas(gcf,[pwd '/figures/FracBlockage1.png'])
saveas(gcf,[pwd '/figures/FracBlockage1'],'epsc')
