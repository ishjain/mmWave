% BlockageNLOSfn.m
% June10: this function implements NLOS blockage

function output = BlockageSimFnNLOS(s_mobility,BS_input,RefInput)
% Update Jun6:BlockageSimFn.m renamed the function and few other parameter's name.
% Update Apr13: Move s_mobility to Simulation.m
% Update Feb 28: Save whole dataBS for all BS and BL density for all aID
%               : Then write separate code for analysis and theoretical
%               plots
% BlockageSim_Feb17.m
% Random Way-point mobility model for blockers
% Simulating Blockage of nT number of BSs.
% Generate time sequence of blocked/unblocked periods
% AND operator on all those nT time sequences.
% FInally, report the blockage freq and duration of all nT BSs
%%frequency = (per sec)count transitions from 0 to 1 devided by simTime
%%duration: (sec) count total # of 1's multiplied by tstep/blockageCount


%----Play-with-values-here--------------------------------------
wannaplot = BS_input.WANNAPLOT; %1;
nB = BS_input.NUM_BL; %number of blokers
nTorig = BS_input.Original_NUM_AP; %Original APs without considering self blockage
rT =BS_input.LOC_AP_DISTANCE; %location of APs
thetaTorig = BS_input.LOC_AP_ANGLE;%location of APs

frac = BS_input.FRACTION;
omega = BS_input.SELF_BL_ANGLE_OMEGA;  

%%Implementing self-blockage
tempInd =  find(thetaTorig>=omega); %These BSs are not blocked by self-blockage
xT = rT(tempInd).*cos(thetaTorig(tempInd));%location of APs (distance)
yT = rT(tempInd).*sin(thetaTorig(tempInd));%location of APs (angle)
nT = length(tempInd); % number of BS not blocked by self-blockage
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
thetaT = thetaTorig(tempInd); %angle from x-axis for BS not blocked by self-bl
simTime = BS_input.SIMULATION_TIME; %sec Total Simulation time
tstep = BS_input.TIME_STEP; %(sec) time step
mu = BS_input.MU; %Expected bloc dur =1/mu

%---------------I am moving this to Sim...m--------
% locT = BS_input.T_EFF_LOCATION; %AP location
% alphaT = BS_input.T_ANGLE;

% s_mobility = Generate_Mobility(s_input);  

%------------------------------------------------------


%---------------NLOS-----------------------------------

nL=RefInput.NUM_REFLECTOR;
ell=RefInput.LOC_REF_DISTANCE;
phi=RefInput.LOC_REF_ANGLE;
d=RefInput.LENGTH_REFLECTOR;
psi=RefInput.ORIENTATION_REF;
R_NLOS= RefInput.RADIUS_COVERAGE_NLOS;

%find blockage prob from theoretical expn for a particular BS and reflector
    
%End point coordinates:
p1 =  rT.*cos(thetaT);
q1 =  rT.*sin(thetaT);
p2 =  ell.*cos(phi);
q2 = ell.*sin(phi);
a =  2*ell.*sin(psi-phi).*sin(psi);
b =  (-2)*ell.*sin(psi-phi).*cos(psi);
OA = sqrt((d./2+ell.*cos(psi-phi)).^2+(ell.*sin(psi-phi)).^2);
OB = sqrt((d./2-ell.*cos(psi-phi)).^2+(ell.*sin(psi-phi)).^2);
angleOAS = psi + atan((d./2+ell.*cos(psi-phi))/(ell.*sin(psi-phi)));
angleOBS = psi - atan((d./2-ell.*cos(psi-phi))/(ell.*sin(psi-phi)));
x1 = OA.*cos(angleOAS);
y1 = OA.*sin(angleOAS);
x2 = OB.*cos(angleOBS);
y2 = OB.*sin(angleOBS);
for iL = 1:nL
for iT = 1:nTorig


%NLOS coverage indicator random variable
d_NLOS(iL,iT) = sqrt(r(iT).^2+4*ell(iL).^2*sin(psi(iL)-phi(iL)).^2+...
    4*r(iT).*ell(iL).*sin(psi(iL)-phi(iL))*sin(theta(iT)-psi(iL)));
d_NLOSapprox = r(iT)+2*ell(iL);

%Indicates whether NLOS path exists (1) or not (0)
IndicatorRV=  ((q1(r,theta)-b(ell,phi,psi)-...
    ((y1(l,phi,d,psi)-b(l,phi,psi))*(p1(r,theta)-a(l,phi,psi)))/(x1(l,phi,d,psi)-a(l,phi,psi))) >=0) && ...
    ((q2(l,phi)-b(l,phi,psi)-...
    ((y2(l,phi,d,psi)-b(l,phi,psi))*(p2(l,phi)-a(l,phi,psi)))/(x2(l,phi,d,psi)-a(l,phi,psi))) >=0) &&...
    (d_NLOS(r,theta,l,phi,psi)<=R_NLOS);


end
end
%----------------------------------------------------

dataBS = cell(nT,1); 
%dataBS contain array of timestamps of blocker arrival for all BSs,
%independent of which blocker is blocking

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
        for indT = 1:nT %for every BS around the UE (outside self-bl zone)
            %The find_blockage_distance() function is written by Ish Jain
            distance_travelled = find_blockage_distance([loc0,loc1],locT(:,indT),thetaT(indT));
            timeToBl = distance_travelled/velocity; %time to blocking event
            timestampBl = start_time+timeToBl; %timestamp of blockage event
            if(distance_travelled>=0 && timestampBl<=simTime)
                %                 data{indB,indT} = [data{indB,indT},start_time+blockage_time];
                dataBS{indT} = [dataBS{indT}, timestampBl];
                
            end
        end
    end
    
end


totaltime = (simTime)/tstep;
binary_seq = zeros(nT,totaltime); %time sequence for every BS
allBl = ones(1,totaltime); %binary seq of all blocked
Tval = tstep:tstep:totaltime*tstep; %run simulation till tdur with step of tstep
if(wannaplot),figure; hold on; end

for indT = 1:nT
    len =length(dataBS{indT});
    dataBS{indT}(2,:) =  exprnd(1/mu,1,len);
end

for indT = 1:nT
    %     blDur  = exprnd(1/mu);
    for timestamp = 1:size(dataBS{indT},2)
        blDur  = ceil(dataBS{indT}(2,timestamp)/tstep);
        blTime = ceil(dataBS{indT}(1,timestamp)/tstep);
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

%%Return now
output=[avgFreq,avgDur,probAllBl,nTorig,nT];

end

