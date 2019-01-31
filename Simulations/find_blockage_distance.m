function dist = find_blockage_distance(locB,locT,alpha)
% Written by Ish Jain
% NYU Tandon School of Engineering
% Date: June 2018

%Function to find the distance of blockage for a particular AP and
%blocker

%Input: blocker location: locB
%locB = |x_prev x_new|
%       |y_prev y_new|
%Input: AP location: locT
%locT = |x|
%       |y|
%Output: dist= distance travelled by blocker from its starting point to the
%point of intersection. In no ointersection, then dist = -1.

%test---------------------
% locT = [2;2];
% locB = [1,1;0,3];
% alpha = pi/4;
wannaplot=0;
%------------------------
rotMat = [cos(alpha),sin(alpha);-sin(alpha),cos(alpha)]; %Rotation Matrix
locBnew = rotMat*locB;
locTnew = rotMat*locT;
xT = locTnew(1);
yT = locTnew(2);
x0B = locBnew(1,1);
y0B = locBnew(2,1);
x1B = locBnew(1,2);
y1B = locBnew(2,2);
frac = y0B/y1B;
xstar = (x0B+frac*x1B)/(1+frac);
if(xstar>=0 && xstar<=xT && y0B*y1B<=0)
   dist =  sqrt((xstar-x0B)^2+y0B^2);
else
    dist = -1;
end

if(wannaplot)
    close all; figure; hold on;
    plot([0,xT],[0,yT],'-r')
    plot([x0B,x1B],[y0B,y1B])
end