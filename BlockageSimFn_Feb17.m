function output = BlockageSimFn_Feb17(s_input,AP_input)
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

R=AP_input.RADIUS_AROUND_UE; %m Radius
nB = s_input.NB_NODES;%4*R^2*rho_b;%=4000; %number of blokers
% nT = AP_input.nT; %number of APs
rhoB=s_input.DENSITY_BL;
V = 1; %m/s Attension!!!!!!!!!!!!!!!!!!!
rhoT = AP_input.DENSITY_AP;
nT = poissrnd(rhoT*pi*R^2)
% nT=1
if(nT==0), output=[0,0,0,0,0,0]; return; end % Dealing zero APs
frac = AP_input.FRACTION;
rT = R*sqrt(rand(nT,1)); %location of APs
alphaT = 2*pi*rand(nT,1);%location of APs
xT = rT.*cos(alphaT);%location of APs
yT = rT.*sin(alphaT);%location of APs
xTfrac = frac*xT; %blockage zone around UE for each APs
yTfrac = frac*yT;
locT = [xTfrac';yTfrac']; %2 rows for x and y, nT columns



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


totaltime = (simTime)/tstep
binary_seq = zeros(nT,totaltime);
allBl = ones(1,totaltime); %binary seq of all blocked
Tval = tstep:tstep:totaltime*tstep; %run simulation till tdur with step of tstep
if(wannaplot),figure; hold on; end

for indT = 1:nT
    %     blDur  = exprnd(1/mu);
    for timestamp = 1:length(dataAP{indT})
        blDur  = ceil(exprnd(1/mu)/tstep);
        blTime = ceil(dataAP{indT}(timestamp)/tstep);
        if(blTime+blDur<=simTime/tstep)%avoid excess duration
            binary_seq(indT, blTime:1:(blTime+blDur))=1;
        end
    end
    allBl = allBl & binary_seq(indT,:);
    if(wannaplot)
        subplot(nT+1,1,indT)
        plot(binary_seq(indT,:), 'lineWidth',4)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
    end
end
if(wannaplot)
    subplot(nT+1,1,nT+1)
    plot(Tval,allBl,'r-', 'lineWidth',4)
    xlabel('Time (sec)')
end
% plot(binary_seq(1,:), 'lineWidth',4)


%%Evaluate frequency and average duration of blockage
avgFreq = sum(diff(allBl)>0)/simTime;
avgDur = sum(allBl)*tstep/sum(diff(allBl)>0);
probAllBl = sum(allBl)*tstep/simTime;

%%Get theoretical values
temp = 2/pi*rhoB*V*frac;
c = temp/mu;

%%This is average frequency of blockage theoretical value
% th_freqBl = 2/pi*rhoB*V*frac*sum(rT); % Theoretical rate of blocking 1 AP
        
a = 1-2*mu./(R*temp) + 2*mu^2./(R^2*temp.^2).*log(1+temp.*R/mu);
th_freqBl = mu*a.*rhoT*pi*R^2.*exp((a-1).*rhoT*pi*R^2);


th_durBl = 1/(nT*mu);
th_probAllBl = exp(-2*pi.*R.*rhoT/c).*(1+c*R).^(2*pi.*rhoT/c^2);

%%Return now
output=[avgFreq,avgDur,probAllBl,th_freqBl,th_durBl,th_probAllBl];


end