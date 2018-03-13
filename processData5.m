%Process Data5 and save the output to simData.csv
%Mar13: processData5.m: Got Data5 from Rajeev 10,000 data but many missing. Use
%updated theoretical result (from Theory.m) and update sim result by also
%finding the conditional probabililities and expectations of freq and dur.
%processData2: for mu=2 and good files

close all
clear
wannaplot=1;
nFiles = 10000;

densityBL = [0.01,0.1,0.2,0.5,0.65];
densityAP = (1:1:10)/10^4;
nBL = length(densityBL);
nAP = length(densityAP);
tempInd = 0;
tempInd2=0;
num0BS = zeros(1,length(densityAP)*length(densityBL));
for i=1:nFiles
    if (exist(strcat('Data5/output',int2str(i),'.csv'))==0)
        continue;
    else
        tempInd=tempInd+1;
        data(:,:,tempInd)=csvread(strcat('Data5/output',int2str(i),'.csv'));
        colNum = find(~any(data(:,:,tempInd),1));
        if(isempty(colNum))
            colNum;
            tempInd2=tempInd2+1;
            data_nn0(:,:,tempInd2) = data(:,:,tempInd);
        end
        num0BS(colNum) = num0BS(colNum)+1;
             
    end
end
num0BS = reshape(num0BS, nBL,nAP);
% newData = reshape(nBL,nAP);

data(isnan(data))=0;
meanVal = reshape(mean(data,3),6,nBL,nAP);

data_nn0(isnan(data_nn0))=0;
meanVal0 = reshape(mean(data_nn0,3),6,nBL,nAP);
%The 6 columns are [avgFreq,avgDur,probAllBl,th_freqBl,th_durBl,th_probAllBl];
avgFreq = squeeze(meanVal(1,:,:));
avgDur= squeeze(meanVal(2,:,:));
probAllBl= squeeze(meanVal(3,:,:));

pB=squeeze(meanVal(3,:,:));
pBgiven=squeeze(meanVal0(3,:,:));
freq=squeeze(meanVal(1,:,:));
freqCond=squeeze(meanVal0(1,:,:));
durCond=squeeze(meanVal0(2,:,:));

allData = [pB,pBgiven,freq,freqCond,durCond];
% csvwrite('simData.csv',allData);
csvwrite('figures2/sim_pB.csv',[densityAP;pB]');
csvwrite('figures2/sim_pBCond.csv',[densityAP;pBgiven]');
csvwrite('figures2/sim_freq.csv',[densityAP;freq]');
csvwrite('figures2/sim_freqCond.csv',[densityAP;freqCond]');
csvwrite('figures2/sim_durCond.csv',[densityAP;durCond]');


if(wannaplot)
    figure(1);
    semilogy(densityAP,pB); 
    ylim([1e-4,1]);title('Marginal prob of Blockage')
    
    figure(2);
    semilogy(densityAP,pBgiven); title('Conditional prob of Bl given n!=0')
    ylim([1e-4,1])
    
    figure(3);
    semilogy(densityAP,freq)
    title('Expected Freq of blockage')
%     ylim([1e-4,1])
        figure(4);
    semilogy(densityAP,freqCond);
    title('Conditional expectation of freq of bl given n!=0')
%     ylim([1e-4,1])
        figure(5);
   semilogy(densityAP,durCond)
    title('Conditional expectation of duration of bl given n!=0')
    
    
end





