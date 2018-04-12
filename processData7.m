%Process Data5 and save the output to figures2/simData csv files
%Mar13: processData5.m: Got Data5 from Rajeev 10,000 data but many missing. Use
%updated theoretical result (from Theory.m) and update sim result by also
%finding the conditional probabililities and expectations of freq and dur.
%processData2: for mu=2 and good files

close all
clear
wannaplot=1;
nFiles = 10000;%3768;

% densityBL = [0.01,0.1,0.2,0.5,0.65];
% densityAP = (1:1:10)/10^4;
densityBL = [0.1,0.2];
densityAP = [10,20,50,100,200,500]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, pi/3, pi/2, 2*pi/3];
mu=2;

nBL = length(densityBL);
nAP = length(densityAP);
nO = length(omegaVal);
tempInd = 0;
% tempInd2=0;
num0BS = zeros(1,nBL*nAP*nO);
num0BS_debug = zeros(1,nBL*nAP*nO);
Directory = 'Data7\';
% Directory = 'rajeevNew\';
for i=1:nFiles
    if (exist(strcat(Directory,'output',int2str(i),'.csv'))==0)
        continue;   
   
    else
        tempInd=tempInd+1;
        data(:,:,tempInd)=csvread(strcat(Directory,'output',int2str(i),'.csv'));
        colNum0 = find(data(5,:,tempInd)==0); %column index of data where num AP =0
        colNum0debug = find(data(4,:,tempInd)==0);
        num0BS_debug(colNum0debug) = num0BS_debug(colNum0debug)+1;
        %         colNum_not0(1,:,tempInd) = find(data(5,:,tempInd)~=0);
        NaNarray = isnan(data(2,:,tempInd));
        %         countNaN(NaNarray) = countNaN(NaNarray)+1;
        %when dur in NaN, replace it by 1/(n\mu);
        data(2,NaNarray,tempInd) = 1./(mu*data(5,NaNarray,tempInd));
        %         colNum2 = find(~any(data(:,:,tempInd),1));
        %         if(isempty(colNum))
        %             colNum;
        %             tempInd2=tempInd2+1;
        %             data_nn0(:,:,tempInd2) = data(:,:,tempInd);
        %         end
        data_nn0(:,:,tempInd) = data(:,:,tempInd);
        num0BS(colNum0) = num0BS(colNum0)+1;
        data(3,colNum0,tempInd) = 1; %when n=0, prob of blockage=1;
    end
end
ratio = num0BS_debug./num0BS;

% num0BS = reshape(num0BS, nBL,nAP,nO);
% newData = reshape(nBL,nAP);

% data(isnan(data))=0;
% countNaNmatrix = isnan(data);
% countNaN = sum(countNaNmatrix,3);
%
%
% %
% meanVal0forDur = sum(data_nn0,3)./(size(data_nn0,3)-countNaN);
% temp = reshape(meanVal0forDur,size(data_nn0,1),nBL,nAP);
% durCond=squeeze(temp(2,:,:));


meanData = mean(data,3);
stdData = std(data,0,3);

meanVal = reshape(meanData,size(data,1),nBL,nAP,nO);
stdVal = reshape(stdData,size(data,1),nBL,nAP,nO);

meanData0 = sum(data_nn0,3)./(size(data_nn0,3)-meshgrid(num0BS,1:size(data,1)));
meanVal0 = reshape(meanData0,size(data,1),nBL,nAP,nO);
% meanVal0 = reshape(mean(data_nn0,3),size(data,1),nBL,nAP);
%The 6 columns are [avgFreq,avgDur,probAllBl,th_freqBl,th_durBl,th_probAllBl];
% avgFreq = squeeze(meanVal(1,:,:,:));
% avgDur= squeeze(meanVal(2,:,:,:));
% pB= squeeze(meanVal(3,:,:,:));

pB=squeeze(meanVal(3,:,:,:));
std_pB = squeeze(stdVal(3,:,:,:));

pBCond=squeeze(meanVal0(3,:,:,:));


freq=squeeze(meanVal(1,:,:,:));
std_freq=squeeze(stdVal(1,:,:,:));

freqCond=squeeze(meanVal0(1,:,:,:));
durCond=squeeze(meanVal0(2,:,:,:));

%%Make it 2D arrays
pB=cat(1,pB(:,:,1),pB(:,:,2))';
std_pB=cat(1,std_pB(:,:,1),std_pB(:,:,2))';

pBCond=cat(1,pBCond(:,:,1),pBCond(:,:,2))';
freq=cat(1,freq(:,:,1),freq(:,:,2))';
freqCond=cat(1,freqCond(:,:,1),freqCond(:,:,2))';
durCond=cat(1,durCond(:,:,1),durCond(:,:,2))';
% pB=cat(1,pB(:,:,1),pB(:,:,2))';
% pBCond=cat(1,pBCond(:,:,1),pBCond(:,:,2))';
% freq=cat(1,freq(:,:,1),freq(:,:,2))';
% freqCond=cat(1,freqCond(:,:,1),freqCond(:,:,2))';
% durCond=cat(1,durCond(:,:,1),durCond(:,:,2))';
%6 rows represents 6 values of AP density
%8 columns has the following legend
% legendArray = {'\lambda_B=0.1,\omega=0','\lambda_B=0.2,\omega=0',...
%     '\lambda_B=0.1,\omega=\pi/3','\lambda_B=0.2,\omega=\pi/3',...
%     '\lambda_B=0.1,\omega=\pi/2','\lambda_B=0.2,\omega=\pi/2',...
%     '\lambda_B=0.1,\omega=2\pi/3','\lambda_B=0.2,\omega=2\pi/3'};
legendArray = {'\lambda_B=0.1,\omega=0','\lambda_B=0.2,\omega=0',...
    '\lambda_B=0.1,\omega=2\pi/3','\lambda_B=0.2,\omega=2\pi/3'};

colTitle= {'lamT','lamB0.1omega0','lamB0.2omega0',...
    'lamB0.1omega2pi3','lamB0.2omega2pi3'};
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,pB])]),...
    'figures2/sim_pB.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,pBCond])]),...
    'figures2/sim_pBCond.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,freq])]),...
    'figures2/sim_freq.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,freqCond])]),...
    'figures2/sim_freqCond.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,durCond])]),...
    'figures2/sim_durCond.csv','writevariablenames',0)

% csvwrite('figures2/sim_pB.csv',[densityAP'*10^4,pB]);
% csvwrite('figures2/sim_pBCond.csv',[densityAP'*10^4,pBCond]);
% csvwrite('figures2/sim_freq.csv',[densityAP'*10^4,freq]);
% csvwrite('figures2/sim_freqCond.csv',[densityAP'*10^4,freqCond]);
% csvwrite('figures2/sim_durCond.csv',[densityAP'*10^4,durCond]);


if(wannaplot)
    figure(1);
    errorbar(meshgrid(densityAP,1:4)',pB,pB-std_pB,pB+std_pB,'o');
    set(gca,'YScale','log');
    ylim([1e-4,1]);title('Marginal prob of Blockage')
    legend(legendArray);
    figure(2);
    semilogy(densityAP,pBCond); title('Conditional prob of Bl given n!=0')
    ylim([1e-4,1])
    legend(legendArray);
    figure(3);
    semilogy(densityAP,freq)
    title('Expected Freq of blockage')
    legend(legendArray);
    ylim([1e-4,1])
    figure(4);
    semilogy(densityAP,freqCond);
    title('Conditional expectation of freq of bl given n!=0')
    ylim([1e-4,1])
    legend(legendArray);
    figure(5);
    semilogy(densityAP,durCond)
    title('Conditional expectation of duration of bl given n!=0')
    legend(legendArray);
    
end





