% Jun6: Renamed to SinumationLOS.m and added some comments
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
if(isempty(aID))
  warning('aID is empty. Replacing it with 1.')  
  aID = '1'; %Runs only for first value of AP density when aID=1
end
rng('shuffle');
wannaplot=0; %Don't plot if running for many loops (else too many plots).
V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
simTime = 3*60*60; %sec Total Simulation time
% Note!!! simTime must be >100ms else the code won't work :)
tstep = 0.0001; %(sec) time step
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
densityBL = [0.005,0.01];
densityAP = [50,100,200,300,400,500]*10^(-6);%(1:1:10)/10^4;
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

% Generate_Mobility function is Copyright (c) 2011, Mathieu Boutin
s_mobility{indB} = Generate_Mobility(s_input{indB});
end
finaldata = zeros(5,length(densityAP),length(densityBL),length(omegaVal));

for indT = 1:length(densityAP)
        rhoT = densityAP(indT);
        nTorig = poissrnd(rhoT*pi*R^2); %original AP number (without self-block)
        rT = R*sqrt(rand(nTorig,1)); %location of APs (distance from origin)
        alphaT = 2*pi*rand(nTorig,1);%location of APs (angle from x-axis)
    for indB = 1:length(densityBL) %for all blockers
        for indO = 1:length(omegaVal) %for all self-blockage angle
            omega = omegaVal(indO);
            rhoB = densityBL(indB);%0.65;%Rajeev calculated central park
            nB = 4*R^2*rhoB;%=4000; %number of blokers
            
            BS_input = struct('WANNAPLOT',wannaplot,...
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
            
            %BlockageSimFn function is written by Ish Jain
            output = BlockageSimFn(s_mobility{indB},BS_input);
            finaldata(:,indT,indB,indO) = output;
            %         output is [avgFreq,avgDur,probAllBl,th_freqBl,th_durBl,th_probAllBl];
        end
    end
end
%Use the code processData9.m to analyze and plot the results
csvwrite(strcat('output',num2str(aID),'.csv'),finaldata)
