%find_blockage_distance.m
%Feb19
%Function to find the distance of blockage for a particular AP and
%blocker

%Input: blocker location: locB
%locB = |x_prev x_new|
%       |y_prev y_new|
%Input: AP location: locT
%locT = |0 x|
%       |0 y|
%Output: dist= distance travelled by blocker from its starting point to the
%point of intersection. In no ointersection, then dist = -1.

%test---------------------
locT = [0,2;0,2];
locB = [1,1;0,3];
wannaplot=1;
%------------------------
xT = locT(1,2);
yT = locT(2,2);
x0B = locB(1,1);
y0B = locB(2,1);
x1B = locB(1,2);
y1B = locB(2,2);
angleT = 45;

if(wannaplot)
    close all; figure; hold on;
    plot([0,xT],[0,yT],'-r')
    plot([x0B,x1B],[y0B,y1B])
end