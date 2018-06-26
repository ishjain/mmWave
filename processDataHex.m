% processDataHex1.m
% June5: Added self-blockage for 60% worst case with 10 BS blocked
% May24: Data with hexagonal cell numerical integration

close all
clear
wannaplot=1;

densityBL = [0.01,0.1];
densityAP = [50,100,150,175,200,225,250,300,350,400,450,500]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, pi/3]; %It will work only for [0,pi/3], so don't change :P

d = sqrt(2./(3*sqrt(3)*densityAP));
Directory = 'DataHex3\';
nBL = length(densityBL);
nAP = length(densityAP);
nOmega = length(omegaVal);
tempi=0;
for i=1:nAP
    
    if (exist(strcat(Directory,'output',int2str(i),'.csv'))==0)
        continue;
        
    else
        tempi=tempi+1;
        data0=csvread(strcat(Directory,'output',int2str(i),'.csv'));
        hex_data_plot(tempi,:) = [densityAP(i) reshape(data0,1,4)];
        hex_data_plot2(tempi,:) = [ d(i) reshape(data0,1,4)];
    end
end

legendArray= {'\lambda_B=0.01,\omega=0^o','\lambda_B=0.1,\omega=0^o',...
    '\lambda_B=0.01,\omega=60^o','\lambda_B=0.1,\omega=60^o',};
colTitle= {'densityAP','lamB0.01omega0','lamB0.1omega0','lamB0.01omega60','lamB0.1omega60'};
writetable(cell2table([colTitle; num2cell(hex_data_plot)]),...
    'figuresHex/pB_Hex.csv','writevariablenames',0)

if(wannaplot);
    h=figure(1);
    semilogy(hex_data_plot(:,1),hex_data_plot(:,2:end), 'LineWidth',2, 'Marker','o');
    hold on;
    plot([0.5,5]*1e-4, [1e-5,1e-5],'k--')
    title('Blockage probability for Hexagonal cell case');
    
    g=legend(legendArray);
    ylim([1e-7,1])
    set(g,'fontsize',13)
    xlabel('AP Density \lambda_T (x 100/km^2)', 'fontsize',13)
    ylabel('Blockage Probability','fontsize',13)
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,[pwd '/figuresHex/pB_Hex.pdf'],'-dpdf','-r0')
    print(h,[pwd '/figuresHex/pB_Hex.png'],'-dpng','-r0')
end



