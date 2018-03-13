%processData4: Use mat file instead of CSV files to debug etc. 
%Got initial 800 and a total of 10000 iterations result.
% Mar 6th: In each iteration, APs selected according to PPP and blockers
% move acc to random way point model. We save the freq and dur of blockage
% of all AP for all blocker as time sequences. We have 10*5=50 such
% sequences for 10 diff AP densities and 5 diff Blocker densities. Now find
% out the different statistics such as blockage of some or all APs, debug
% etc

clear
close all;
wannaplot=0;
% V = 1; %velocity m/s
% hb = 1.8;
% hr = 1.4;
% ht = 6;
% frac = (hb-hr)/(ht-hr);
% mu = 2; %Expected bloc dur =1/mu
% R = 100; %m Radius
simTime = 60*10; %sec Total Simulation time
tstep = 0.0001; %(sec) time step

densityBL = [0.01,0.1,0.2,0.5,0.65];
densityAP = (1:1:10)/10^4;
aIDval = 10:1:100;
naID = length(aIDval);

count0BS = zeros(length(densityAP),length(densityBL));
freq = zeros(length(densityAP),length(densityBL));
dur = zeros(length(densityAP),length(densityBL));
probAllBl = zeros(length(densityAP),length(densityBL));
for indT = 1:length(densityAP)
    for indB = 1:length(densityBL)
        
        rhoB = densityBL(indB);
        rhoT = densityAP(indT);
        for aIDind=1:length(aIDval);
            aID = aIDval(aIDind);
            filename = strcat('dataAP_',num2str(aID),...
                '_',num2str(indB),...
                '_',num2str(indT),'.mat');
            
            %load dataAP which is structure with nAP elements
            %each element is 2D array, first array for freq (0 to 600sec)
            %and second for duration (1/mu = 1/2sec) timestamp.
            load(fullfile('thesis',filename))
            nT = size(dataAP,1);
            
            %%Here my goal is to deal witj nT = 0;
            if(nT==0)
                count0BS(indT,indB) = count0BS(indT,indB)+1;
%                 freq(indT,indB) = freq(indT,indB) + 0;
                %dur(indT,indB) = dur(indT,indB) + sum(allBl)*tstep/sum(diff(allBl)>0);
%                 probAllBl(indT,indB) = probAllBl(indT,indB) +1;
                
            else
            
            %My purpose is to debug if there is overlapping of two blockage
            %event for single AP.
%             startTime=sort(dataAP{1}(1,:));
%             
%             endTime = startTime+dataAP{1}(2,:);
%             shiftStartTime = startTime(2:end);
%             shiftEndTime = endTime(1:end-1);
%             indx = find((shiftStartTime-shiftEndTime)<0);
            
            %%CopyPaste code from BlockageSimFn_Feb17.m
            
            
            totaltime = (simTime)/tstep;
            binary_seq = zeros(nT,totaltime);
            allBl = ones(1,totaltime); %binary seq of all blocked
            Tval = tstep:tstep:totaltime*tstep; %run simulation till tdur with step of tstep
            
            if(wannaplot),figure; hold on; 
            end
            for n = 1:nT
                %     blDur  = exprnd(1/mu);
                for timestamp = 1:length(dataAP{n})
                    blDur  = ceil(dataAP{n}(2,timestamp)/tstep);
                    blTime = ceil(dataAP{n}(1,timestamp)/tstep);
                    if(blTime+blDur<=simTime/tstep)%avoid excess duration
                        binary_seq(n, blTime+1:1:(blTime+blDur))=...
                            binary_seq(n, blTime+1:1:(blTime+blDur))+1;
                    end
                end
                allBl = allBl & binary_seq(n,:);
                if(wannaplot)
                    range=1:600000;
                    subplot(nT+1,1,n)
                    plot(binary_seq(n,range), 'lineWidth',4)
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                    title('Individual blockage of single AP')
                end
            end
            if(wannaplot)
                subplot(nT+1,1,nT+1)
                plot(Tval(range),allBl(range),'r-', 'lineWidth',4)
                xlabel('Time (sec)')
                title('Simultaneous Blockage of all APs')
            end
            % plot(binary_seq(1,:), 'lineWidth',4)
            
            
            %%Evaluate frequency and average duration of blockage
            freq(indT,indB) = freq(indT,indB)+sum(diff(allBl)>0)/simTime;
            dur(indT,indB) = dur(indT,indB) + sum(allBl)*tstep/sum(diff(allBl)>0);
            probAllBl(indT,indB) = probAllBl(indT,indB) + sum(allBl)*tstep/simTime;
            end

            
            %%Return now
%             output=[avgFreq,avgDur,probAllBl,th_freqBl,th_durBl,th_probAllBl];
            
            
        end
        
    end
    
end
%Take average now
pBAvg= (probAllBl+count0BS)/naID;%!!!!wrong
pBgivenAvg=probAllBl./(naID-count0BS);
freqAvg = freq./(naID);
freqGivenAvg = freq./(naID-count0BS);
durAvg = dur./(naID-count0BS);




