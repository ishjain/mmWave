% Simulation_Feb20
%Apr3: Generate BS number and location once here and run for various
%blocker density and various omega. Earlier it was di=one in BlockageSimFn

% Mar14: Added code for self blockage (omega) and modified bl/AP densities
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
wannaplot=1;
V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 5;
frac = (hb-hr)/(ht-hr);
simTime = 10*60*60; %sec Total Simulation time
tstep = 0.0001; %(sec) time step
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
densityBL = [0.01,0.02];
densityAP = [50,100,200,300,400,500,600]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, pi/3];

s_input = cell(1,2);
s_mobility = cell(1,2);
for indB=1:length(densityBL)
s_input{indB} = struct('V_POSITION_X_INTERVAL',[-R R],...%(m)
    'V_POSITION_Y_INTERVAL',[-R R],...%(m)
    'V_SPEED_INTERVAL',[V V],...%(m/s)
    'V_PAUSE_INTERVAL',[0 0],...%pause time (s)
    'V_WALK_INTERVAL',[1.00 60.00],...%walk time (s)
    'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
    'SIMULATION_TIME',simTime,...%(s)
    'NB_NODES',4*R^2*densityBL(indB));
s_mobility{indB} = Generate_Mobility(s_input{indB});
end
finaldata = zeros(5,length(densityAP),length(densityBL),length(omegaVal));

for indT = 1:length(densityAP)
        rhoT = densityAP(indT);
        nTorig = poissrnd(rhoT*pi*R^2); %original AP number (without self-block)
        rT = R*sqrt(rand(nTorig,1)); %location of APs
        alphaT = 2*pi*rand(nTorig,1);%location of APs
    for indB = 1:length(densityBL)
        
        
        for indO = 1:length(omegaVal)
            omega = omegaVal(indO);
            rhoB = densityBL(indB);%0.65;%Rajeev calculated central park
            nB = 4*R^2*rhoB;%=4000; %number of blokers
            
            AP_input = struct('WANNAPLOT',wannaplot,...
                'RADIUS_AROUND_UE',R,...
                'DENSITY_AP',rhoT,...
                'SIMULATION_TIME',simTime,...
                'TIME_STEP',tstep,...
                'MU',mu,...
                'FRACTION',frac,...
                'SELF_BL_ANGLE_OMEGA',omega,...
                'Original_NUM_AP',nTorig,...
                'LOC_AP_DISTANCE', rT,...
                'LOC_AP_ANGLE',alphaT,...
                'NUM_BL',nB);
            output = BlockageSimFn_Feb17(s_mobility{indB},AP_input);
            finaldata(:,indT,indB,indO) = output;
            %         output is [avgFreq,avgDur,probAllBl,th_freqBl,th_durBl,th_probAllBl];
        end
    end
end

csvwrite(strcat('output',num2str(aID),'.csv'),finaldata)