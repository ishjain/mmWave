% Simulation_Feb20

% Update Feb20: Generate AP location using PPP

% Update Feb17: Transferred to BlockageSimFn_Feb17.m
% Random Way-point mobility model for blockers
% Simulating Blockage of nT number of APs.
% Generate time sequence of blocked/unblocked periods

close all;
clear;

%----Play-with-values---------------------------------------
aID = getenv('SLURM_ARRAY_TASK_ID')
rng('shuffle');
wannaplot=0;
V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 6;
frac = (hb-hr)/(ht-hr);
simTime = 60*10; %sec Total Simulation time
tstep = 0.0001; %(sec) time step
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
densityBL = [0.01,0.1,0.2,0.5,0.65];
densityAP = (1:.5:10)/10^4;



% MAX=200;
% for iter=1:MAX
for indB = 1:length(densityBL)
    for indT = 1:length(densityAP)
        
        rhoB = densityBL(indB);%0.65;%Rajeev calculated central park
        nB = 4*R^2*rhoB;%=4000; %number of blokers
        
        rhoT = densityAP(indT);

        s_input = struct('V_POSITION_X_INTERVAL',[-R R],...%(m)
            'V_POSITION_Y_INTERVAL',[-R R],...%(m)
            'V_SPEED_INTERVAL',[V V],...%(m/s)
            'V_PAUSE_INTERVAL',[0 0],...%pause time (s)
            'V_WALK_INTERVAL',[1.00 60.00],...%walk time (s)
            'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
            'SIMULATION_TIME',simTime,...%(s)
            'NB_NODES',nB,...
            'DENSITY_BL',rhoB);
        
        AP_input = struct('WANNAPLOT',wannaplot,...
            'RADIUS_AROUND_UE',R,...
            'DENSITY_AP',rhoT,...
            'SIMULATION_TIME',simTime,...
            'TIME_STEP',tstep,...
            'MU',mu,...
            'FRACTION',frac);
        output = BlockageSimFn_Feb17(s_input,AP_input);
        finaldata(:,indB,indT) = output;
%         output is [avgFreq,avgDur,probAllBl,th_freqBl,th_durBl,th_probAllBl];
    end
% end
% finaldata(:,:,iter) = output;

end

csvwrite(strcat('output',num2str(aID),'.csv'),finaldata)