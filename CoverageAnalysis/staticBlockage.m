% Written by Ish Jain
% NYU Tandon School of Engineering
% Date: June 2018
% 
% Description:
% Finding coverage probability considering static and Self blockage
% Coverage is as there is atleast one BS not blocked by Self or static blockage.

clear
close all
wannaplot=1;

R = 100; %m Radius
densityD = [1e-8,0.001]; %D = static blockage
densityBS = (50:500)*10^-6;%[50,100,200,300,400,500]*10^(-6);%BS = basestation
omegaVal = [0, pi/3]; %omega values

El = 10; %m
Ew = 10; %m

for iT = 1:length(densityBS)
    tempind=0;
    for iD = 1:length(densityD)
        for iO = 1:length(omegaVal)
            tempind=tempind+1;
            lamT = densityBS(iT);
            lamD = densityD(iD);
            omega = omegaVal(iO);
            beta(iD) = 2/pi*lamD*(El+Ew);
            beta0(iD) = lamD*El*Ew;
            p(iO) = 1- omega/2/pi;
            b(iD) = 2*exp(-beta0(iD))/(beta(iD)^2*R^2) * ...
                (1-(1+beta(iD)*R)*exp(-beta(iD)*R));
            %coverage probability given by static and self blockage
            Pcoverage(iT,tempind) = 1-exp(-b(iD)*p(iO)*lamT*pi*R^2);
            
            
            colTitle{1}='lamT';
            
            
            colTitle{tempind+1} = strcat(...
                'lamD',num2str(lamD*1e4),'omega',num2str(omega*360/2/pi));
            legendArray{tempind} = strcat(...
                'lamD',num2str(lamD*1e4),'omega',num2str(omega*360/2/pi));
        end
    end
end
colorVal = {'c', 'm', 'r', 'b'};
figure;hold on;
tempind=0;
% legendArray = {'\lambda_d=0, \omega=0','\lambda_d=0, \omega = 90',...
%     '\lambda_d=1e-3, \omega=0','\lambda_d=1e-3, \omega = 90'};
for iD = 1:length(densityD)
    for iO=1:length(omegaVal)
        tempind=tempind+1;
        plot(densityBS,Pcoverage(:,tempind),'color',colorVal{tempind});
    end
end
% colTitle = {'lamT','\lamD0omega0','lamD0omega90',...
%     'lambd=1e-3, \omega=0','\lambda_d=1e-3, \omega = 90' }
writetable(cell2table([colTitle; num2cell([densityBS'*10^4, Pcoverage])]),...
    'figures/coverageLOS.csv','writevariablenames',0)

g=legend(legendArray);
set(g,'fontsize',13)
xlabel('AP Density \lambda_T (x 100/km^2)', 'fontsize',13)
ylabel('P(Coverage)','fontsize',13)