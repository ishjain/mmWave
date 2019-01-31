function test_Animate(s_mobility,s_input,time_step)
% Copyright (c) 2011, Mathieu Boutin
% All rights reserved.
% Modified by Ish Jain
% NYU Tandon School of Engineering
% Date: June 2018
v_t = 0:time_step:s_input.SIMULATION_TIME;

for nodeIndex = 1:s_mobility.NB_NODES
    %Simple interpolation (linear) to get the position, anytime.
    %Remember that "interp1" is the matlab function to use in order to
    %get nodes' position at any continuous time.
    vs_node(nodeIndex).v_x = interp1(s_mobility.VS_NODE(nodeIndex).V_TIME,s_mobility.VS_NODE(nodeIndex).V_POSITION_X,v_t);
    vs_node(nodeIndex).v_y = interp1(s_mobility.VS_NODE(nodeIndex).V_TIME,s_mobility.VS_NODE(nodeIndex).V_POSITION_Y,v_t);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;
%Ish add the BSs and links
hb = 1.8;
hr = 1.4;
ht = 6;
R=100;
frac = (hb-hr)/(ht-hr);
lamT = 2e-4; %%%!!!!!!!!!!change here!!!!!!!!!!!!!!!
nB = s_input.NB_NODES;
nT = poissrnd(lamT*pi*R^2);
rT = R*sqrt(rand(nT,1)); %location of APs
alphaT = 2*pi*rand(nT,1);%location of APs
xT = rT.*cos(alphaT);%location of APs
yT = rT.*sin(alphaT);%location of APs
xTfrac = frac*xT; %blockage zone around UE for each APs
yTfrac = frac*yT;
scatter(xT,yT,30,'r^')

plot([zeros(length(xT),1),xT]',[zeros(length(yT),1),yT]','g')
plot([zeros(length(xTfrac),1),xTfrac]',[zeros(length(yTfrac),1),yTfrac]',...
    'r','linewidth',1)
th=0:0.01:2*pi; xx = R*cos(th); yy = R*sin(th);
plot(xx,yy)
scatter(0,0,20,'s')

for nodeIndex = 1:s_mobility.NB_NODES
    vh_node_pos(nodeIndex) = plot(vs_node(nodeIndex).v_x(1),...
        vs_node(nodeIndex).v_y(1),'.','color',[0.3 0.3 1],...
        'MarkerSize',10);
end
title(cat(2,'Simulation time (sec): ',num2str(s_mobility.SIMULATION_TIME)));
xlabel('X (meters)');
ylabel('Y (meters)');
xlim([-100,100])
ylim([-100,100])
name=strcat('Simulation: nT=',num2str(nT),'  nB=',...
    num2str(nB), '  hT=',num2str(ht));
title(name);
ht = text(-100,90,cat(2,'Time (sec) = 0'));
%     ht = text(min(vs_node(1).v_x),max(vs_node(1).v_y),cat(2,'Time (sec) = 0'));
%     axis([min(vs_node(1).v_x) max(vs_node(1).v_x) min(vs_node(1).v_y) max(vs_node(1).v_y)]);
hold off;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!change file name
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
v = VideoWriter('simulation1','MPEG-4');
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
v.FrameRate=10;
%     v.CompressionRatio = 3;
for timeIndex = 1:length(v_t);
    t = v_t(timeIndex);
    set(ht,'String',cat(2,'Time (sec) = ',num2str(t,4)));
    for nodeIndex = 1:s_mobility.NB_NODES
        set(vh_node_pos(nodeIndex),'XData',vs_node(nodeIndex).v_x(timeIndex),'YData',vs_node(nodeIndex).v_y(timeIndex));
    end
    set(gcf,'paperunits','inches')
    set(gcf,'position',[0 0 550 500])
    drawnow;
    %Ish save the movie
    
    ax = gca;
    
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    F(timeIndex) = getframe(ax,rect);
    
end
open(v)
writeVideo(v,F);
close(v)
end