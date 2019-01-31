
%Process Simulation Data collected from NYU HPC and save the output to 
% figures2/*csv files as well as plot them for visualizing


close all
clear
wannaplot=1;
nFiles = 10000;%3768;

densityBL = [0.005,0.01];
densityAP = [50,100,200,300,400,500]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, pi/3];
% densityBL = [0.01,0.1,0.2];
% densityAP = [50,100,200,300,400,500,600]*10^(-6);%(1:1:10)/10^4;
% omegaVal = [0, pi/3, pi/2];
mu=2;

nBL = length(densityBL);
nAP = length(densityAP);
nO = length(omegaVal);
tempInd = 0;
% tempInd2=0;
num0BS = zeros(1,nBL*nAP*nO);
num0BS_debug = zeros(1,nBL*nAP*nO);
Directory = {'Data9\Part1\','Data9\Part2\'};
% Directory = 'rajeevNew\';
% curr_pos=zeros(1,nBL*nAP*nO);
datanew = cell(1,nBL*nAP*nO);
debugarr = zeros(1,nBL*nAP*nO);
for dir=1:2
    for i=1:nFiles
        if (exist(strcat(Directory{dir},'output',int2str(i),'.csv'))==0)
            continue;
            
        else
            tempInd=tempInd+1;
            %         if (i<=4000)
            %             data0=csvread(strcat(Directory,'output',int2str(i),'.csv'));
            %
            %             %!!!!!!!!my bad!!!!!!!!!!!!!!
            %             data1 = reshape(data0,[5,7,7,3]);
            %             data2 = data1(:,:,1:3,:);
            %             data(:,:,tempInd) = reshape(data2,[5,63]);
            %             %!!!!!!!!!!corrected!!!!!!!!!!!!!!!
            %         else
            data(:,:,tempInd)=csvread(strcat(Directory{dir},'output',int2str(i),'.csv'));
            
            %         end
            colNum0 = find(data(5,:,tempInd)==0); %column index of data where num AP =0
            %         colNum0debug = find(data(4,:,tempInd)==0);
            %         num0BS_debug(colNum0debug) = num0BS_debug(colNum0debug)+1;
            %         colNum_not0(1,:,tempInd) = find(data(5,:,tempInd)~=0);
            NaNarray = isnan(data(2,:,tempInd));
            debugarr = debugarr+NaNarray;
            %         countNaN(NaNarray) = countNaN(NaNarray)+1;
            %when dur in NaN, replace it by 1/(n\mu);
            
            debugdata(:,:,tempInd) = data(:,:,tempInd);
            data(2,NaNarray,tempInd) = 1./(mu*data(5,NaNarray,tempInd));
            
            data_nn0(:,:,tempInd) = data(:,:,tempInd);
            num0BS(colNum0) = num0BS(colNum0)+1;
            data(3,colNum0,tempInd) = 1; %when n=0, prob of blockage=1;
            
            for ii = 1:size(data,2)
                if(data(5,ii,tempInd)~=0)
                    %                 curr_pos(ii) = curr_pos(ii)+1;
                    datanew{ii} = [datanew{ii},data(1:3,ii,tempInd)];
                end
                
            end
        end
    end
end
for ii = 1:size(data,2)
    
    meanval(:,ii) =  mean(datanew{ii},2);
    stdval(:,ii) =  std(datanew{ii},0,2);
    size_data(ii) = size(datanew{ii},2);
end

confidense_int = 1.96*stdval./meshgrid(sqrt(size_data),1:size(stdval,1));
res = [meanval; confidense_int];
res1 = reshape(res,size(res,1),nAP,nBL,nO);


% meanData = mean(data,3);
% meanVal = reshape(meanData,size(data,1),nAP,nBL,nO);

% meanData0 = sum(data_nn0,3)./( size(data_nn0,3)-meshgrid(num0BS,1:size(data,1)) );
% meanVal0 = reshape(meanData0,size(data,1),nAP,nBL,nO);
% confidense_int_new = reshape(confidense_int,size(data,1),nAP,nBL,nO);


for jj = 1:size(res,1) %3 for freq, dur, pB and 3 for their intrvl
    res2{jj} = squeeze(res1(jj,:,:,:));
    res3{jj} = reshape(res2{jj},nAP,nBL*nO);
    
end
pBCond = [res3{3},res3{6}]; %[mean;intrvl]
freqCond = [res3{1},res3{4}];
durCond = [res3{2},res3{5}]*1000; % convert to ms;
densityAP = [50,100,200,300,400,500]*10^(-6);%(1:1:10)/10^4;

% pB=squeeze(meanVal(3,:,:,:));
% % std_pB = squeeze(stdVal(3,:,:,:));
% freq=squeeze(meanVal(1,:,:,:));
% % std_freq=squeeze(stdVal(1,:,:,:));
% freqCond=squeeze(meanVal0(1,:,:,:));
% durCond=squeeze(meanVal0(2,:,:,:));
% pBCond=squeeze(meanVal0(3,:,:,:));
% freqCond_int = squeeze(confidense_int_new(1,:,:,:));
% durCond_int = squeeze(confidense_int_new(2,:,:,:));
% pBCond_int = squeeze(confidense_int_new(3,:,:,:));
% %%Make it 2D arrays
% pB=reshape(pB(1:end-1,:,[1,2]),length(densityAP),6);
% pBCond=reshape(pBCond(1:end-1,:,[1,2]),length(densityAP),6);
% freq=reshape(freq(1:end-1,:,[1,2]),length(densityAP),6);
% freqCond=reshape(freqCond(1:end-1,:,[1,2]),length(densityAP),6);
% durCond=reshape(durCond(1:end-1,:,[1,2]),length(densityAP),6);





legendArray= {'lamB0.005omega0','lamB0.01omega0',...
    'lamB0.005omega60','lamB0.01omega60'};
colTitle= {'lamT','lamB0.005omega0','lamB0.01omega0',...
    'lamB0.005omega60','lamB0.01omega60',...
    'lamB0.005omega0int','lamB0.01omega0int',...
    'lamB0.005omega60int','lamB0.01omega60int'};
% writetable(cell2table([colTitle; num2cell([densityAP'*10^4,pB])]),...
%     'figures2/sim_pB.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,pBCond])]),...
    'figures2/sim_pBCond2.csv','writevariablenames',0)
% writetable(cell2table([colTitle; num2cell([densityAP'*10^4,freq])]),...
%     'figures2/sim_freq.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,freqCond])]),...
    'figures2/sim_freqCond2.csv','writevariablenames',0)
writetable(cell2table([colTitle; num2cell([densityAP'*10^4,durCond])]),...
    'figures2/sim_durCond2.csv','writevariablenames',0)

% csvwrite('figures2/sim_pB.csv',[densityAP'*10^4,pB]);
% csvwrite('figures2/sim_pBCond.csv',[densityAP'*10^4,pBCond]);
% csvwrite('figures2/sim_freq.csv',[densityAP'*10^4,freq]);
% csvwrite('figures2/sim_freqCond.csv',[densityAP'*10^4,freqCond]);
% csvwrite('figures2/sim_durCond.csv',[densityAP'*10^4,durCond]);


if(wannaplot)
    %     figure(1);
    %     semilogy(densityAP,pB);
    %     set(gca,'YScale','log');
    % %     ylim([1e-4,1]);title('Marginal prob of Blockage')
    %     legend(legendArray);
    figure(2);
    semilogy(densityAP,pBCond(:,1:4)); title('Conditional prob of Bl given n!=0')
    %     ylim([1e-4,1])
    legend(legendArray);
    
    
    %     figure(3);
    %     semilogy(densityAP,freq)
    %     title('Expected Freq of blockage')
    %     legend(legendArray);
    % %     ylim([1e-4,1])
    figure(4);
    semilogy(densityAP,freqCond(:,1:4));
    title('Conditional expectation of freq of bl given n!=0')
    %     ylim([1e-4,1])
    legend(legendArray);
    
    figure(5);
    semilogy(densityAP,durCond(:,1:4))
    title('Conditional expectation of duration of bl given n!=0')
    legend(legendArray);
    
end





