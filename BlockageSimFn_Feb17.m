function output = BlockageSimFn_Feb17(s_input,AP_input)
% Update Feb 28: Save whole DataAP for all AP and BL density for all aID
%               : Then write separate code for analysis and theoretical
%               plots
% BlockageSim_Feb17
% Random Way-point mobility model for blockers
% Simulating Blockage of nT number of APs.
% Generate time sequence of blocked/unblocked periods
% AND operator on all those nT time sequences.
% FInally, report the blockage freq and duration of all nT APs
%%frequency = (per sec)count transitions from 0 to 1 devided by simTime
%%duration: (sec) count total # of 1's multiplied by tstep/blockageCount


%----Play-with-values-here--------------------------------------
wannaplot = AP_input.WANNAPLOT; %1;

nB = s_input.NB_NODES;%4*R^2*rho_b;%=4000; %number of blokers

nTorig = AP_input.Original_NUM_AP;
rT =AP_input.LOC_AP_DISTANCE; %location of APs
alphaTorig = AP_input.LOC_AP_ANGLE;%location of APs

frac = AP_input.FRACTION;
omega = AP_input.SELF_BL_ANGLE_OMEGA;  
% psirand = 
tempInd =  find(alphaTorig>=omega);
xT = rT(tempInd).*cos(alphaTorig(tempInd));%location of APs
yT = rT(tempInd).*sin(alphaTorig(tempInd));%location of APs
nT = length(tempInd);
% nT=0
if(nT==0)
    output=[0,0,0,nTorig,nT];   
    return;
end % Dealing zero APs
   
% xT = rT.*cos(alphaT);%location of APs
% yT = rT.*sin(alphaT);%location of APs
xTfrac = frac*xT; %blockage zone around UE for each APs
yTfrac = frac*yT;
locT = [xTfrac';yTfrac']; %2 rows for x and y, nT columns
alphaT = alphaTorig(tempInd);


simTime = AP_input.SIMULATION_TIME; %sec Total Simulation time
tstep = AP_input.TIME_STEP; %(sec) time step
mu = AP_input.MU; %Expected bloc dur =1/mu
% locT = AP_input.T_EFF_LOCATION; %AP location
% alphaT = AP_input.T_ANGLE;

s_mobility = Generate_Mobility(s_input);

%------------------------------------------------------

dataAP = cell(nT,1); %contain array of timestamps for all APs no matter which blocker

for indB = 1:nB %for every blocker
    for iter =1:(length(s_mobility.VS_NODE(indB).V_POSITION_X)-1)
        % for every time blocker changes direction
        loc0 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter);...
            s_mobility.VS_NODE(indB).V_POSITION_Y(iter)];
        loc1 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter+1);...
            s_mobility.VS_NODE(indB).V_POSITION_Y(iter+1)];
        start_time = s_mobility.VS_NODE(indB).V_TIME(iter);
        velocity = sqrt((s_mobility.VS_NODE(indB).V_SPEED_X(iter))^2+ ...
            (s_mobility.VS_NODE(indB).V_SPEED_Y(iter))^2);
        for indT = 1:nT %for every AP
            distance_travelled = find_blockage_distance([loc0,loc1],locT(:,indT),alphaT(indT));
            timeToBl = distance_travelled/velocity; %time to blocking event
            timestampBl = start_time+timeToBl; %timestamp of blockage event
            if(distance_travelled>=0 && timestampBl<=simTime)
                %                 data{indB,indT} = [data{indB,indT},start_time+blockage_time];
                dataAP{indT} = [dataAP{indT}, timestampBl];
                
            end
        end
    end
    
end


totaltime = (simTime)/tstep;
binary_seq = zeros(nT,totaltime);
allBl = ones(1,totaltime); %binary seq of all blocked
Tval = tstep:tstep:totaltime*tstep; %run simulation till tdur with step of tstep
if(wannaplot),figure; hold on; end

for indT = 1:nT
    len =length(dataAP{indT});
    dataAP{indT}(2,:) =  exprnd(1/mu,1,len);
end

% indT = plot_input.indT;
% indB = plot_input.indB;
% aID = plot_input.aID;
% save(strcat('dataAP_',num2str(aID),...
%     '_',num2str(indB),...
%     '_',num2str(indT),'.mat'),'dataAP')
% csvwrite(strcat('output',num2str(aID),'.csv'),finaldata)
for indT = 1:nT
    %     blDur  = exprnd(1/mu);
    for timestamp = 1:size(dataAP{indT},2)
        blDur  = ceil(dataAP{indT}(2,timestamp)/tstep);
        blTime = ceil(dataAP{indT}(1,timestamp)/tstep);
        if(blTime+blDur<=simTime/tstep)%avoid excess duration
            binary_seq(indT, blTime+1:1:(blTime+blDur))=binary_seq(indT, blTime+1:1:(blTime+blDur))+1;
        end
    end
    allBl = allBl & binary_seq(indT,:);
    if(wannaplot)
        subplot(nT+1,1,indT)
        plot(binary_seq(indT,1:10/tstep), 'lineWidth',4)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
    end
end
if(wannaplot)
    subplot(nT+1,1,nT+1)
    plot(Tval(1:10/tstep),allBl(1:10/tstep),'r-', 'lineWidth',4)
    xlabel('Time (sec)')
end
% plot(binary_seq(1,:), 'lineWidth',4)


%%Evaluate frequency and average duration of blockage
avgFreq = sum(diff(allBl)>0)/simTime;
avgDur = sum(allBl)*tstep/sum(diff(allBl)>0);
probAllBl = sum(allBl)*tstep/simTime;

% %%Get theoretical values
% temp = 2/pi*rhoB*V*frac;
% c = temp/mu;
% 
% %%This is average frequency of blockage theoretical value
% % th_freqBl = 2/pi*rhoB*V*frac*sum(rT); % Theoretical rate of blocking 1 AP
%         
% a = 1-2*mu./(R*temp) + 2*mu^2./(R^2*temp.^2).*log(1+temp.*R/mu);
% th_freqBl = mu*a.*rhoT*pi*R^2.*exp((a-1).*rhoT*pi*R^2);
% th_probAllBl = exp(-(1-a)*rhoT);
% 
% th_durBl = 1/(nT*mu);
% % th_probAllBl = exp(-2*pi.*R.*rhoT/c).*(1+c*R).^(2*pi.*rhoT/c^2);

%%Return now
output=[avgFreq,avgDur,probAllBl,nTorig,nT];


end