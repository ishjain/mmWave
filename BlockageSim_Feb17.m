%BlockageSim_Feb17
% Random Way-point mobility model for blockers
% Simulating Blockage of nT number of APs.
% Generate time sequence of blocked/unblocked periods

close all;
clear;

rho_b = 0.65;%Rajeev calculated central park
V=1;
hb = 1.8;
hr = 1.4;
ht = 6;
frac = (hb-hr)/(ht-hr);
temp = 2/pi*rho_b*V*frac;
R=100; %m Radius
distValues = 1:10:100;
tstep = 0.001; %i.e. 1 ms
tdur = 10; %simulation duration
Tval = 0:tstep:tdur; %run simulation till tdur with step of tstep
nB=  100; %number of blokers
nT = 10; %number of APs
rT = R*sqrt(rand(nT,1)); %location of APs
alphaT = 2*pi*rand(nT,1);%location of APs
xT = rT.*cos(alphaT);%location of APs
yT = rT.*sin(alphaT);%location of APs
data = cell(nB,nT); %contans timestamp of blockage of all APs by all blockers
s_input = struct('V_POSITION_X_INTERVAL',[-100 100],...%(m)
                 'V_POSITION_Y_INTERVAL',[-100 100],...%(m)
                 'V_SPEED_INTERVAL',[1 1],...%(m/s)
                 'V_PAUSE_INTERVAL',[0 0],...%pause time (s)
                 'V_WALK_INTERVAL',[1.00 180.00],...%walk time (s)
                 'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
                 'SIMULATION_TIME',500,...%(s)
                 'NB_NODES',100);
s_mobility = Generate_Mobility(s_input);

for iter =1:1
    for indB = 1:nB %for every blocker
       loc0 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter);...
           s_mobility.VS_NODE(indB).V_POSITION_Y(iter)]; 
       loc1 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter+1);...
           s_mobility.VS_NODE(indB).V_POSITION_Y(iter+1)];
       start_time = s_mobility.VS_NODE(indB).V_TIME(iter);
       velocity = sqrt((s_mobility.VS_NODE(indB).V_SPEED_X(iter))^2+ ...
           (s_mobility.VS_NODE(indB).V_SPEED_Y(iter))^2);
       for indT = 1:nT %for every AP
           locT = [xT(indT);yT(indT)];
           
           distance_travelled = find_blockage_distance([loc0,loc1],locT,alphaT(indT));
           if(distance_travelled>=0)
               blockage_time = distance_travelled/velocity;
               
               data{indB,indT} = 1;
           end
       end
    end
    
end
