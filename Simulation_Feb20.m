% Simulation_Feb20

% Update Feb20: Generate AP location using PPP

% Update Feb17: Transferred to BlockageSimFn_Feb17.m
% Random Way-point mobility model for blockers
% Simulating Blockage of nT number of APs.
% Generate time sequence of blocked/unblocked periods

close all;
clear;

%----Play-with-values---------------------------------------

V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 6;

R = 100; %m Radius
densityBL = [0.01,0.1,0.2,0.5,0.65];
rho_b = 0.01;%0.65;%Rajeev calculated central park
nB = 4*R^2*rho_b;%=4000; %number of blokers


nT = 3; %number of APs
simTime = 60*10; %sec Total Simulation time 
tstep = 0.01; %(sec) time step
mu = 5; %Expected bloc dur =1/mu
%------------------------------------------------------

frac = (hb-hr)/(ht-hr);
% temp = 2/pi*rho_b*V*frac; %for theoretical only
rT = R*sqrt(rand(nT,1)); %location of APs
alphaT = 2*pi*rand(nT,1);%location of APs
xT = rT.*cos(alphaT);%location of APs
yT = rT.*sin(alphaT);%location of APs
xTfrac = frac*xT; %blockage zone around UE for each APs
yTfrac = frac*yT;
locT = [xTfrac';yTfrac']; %2 rows for x and y, nT columns
% data = cell(nB,nT); %contans timestamp of blockage of all APs by all blockers

s_input = struct('V_POSITION_X_INTERVAL',[-R R],...%(m)
    'V_POSITION_Y_INTERVAL',[-R R],...%(m)
    'V_SPEED_INTERVAL',[V V],...%(m/s)
    'V_PAUSE_INTERVAL',[0 0],...%pause time (s)
    'V_WALK_INTERVAL',[1.00 60.00],...%walk time (s)
    'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
    'SIMULATION_TIME',simTime,...%(s)
    'NB_NODES',nB);

iterMax =100;
for iter =1:iterMax

AP_input = struct('RADIUS',R,...
    'T_EFF_LOCATION',locT,...
    'T_ANGLE',alphaT,...
    'nT',nT,...
    'SIMULATION_TIME',simTime,...
    'TIME_STEP',tstep,...
    'MU',mu);
[avgFreq,avgDur] = BlockageSimFn_Feb17(s_input,AP_input);

end