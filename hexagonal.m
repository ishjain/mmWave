%May 20: hexagon
% numerical analysis with hexagonal cells
clear

 syms r theta D;

% D = 10;

ri=[repmat(sqrt(3)*D,1,6) repmat(3*D,1,6) repmat(2*sqrt(3)*D,1,6)];
phi = [pi/6, -pi/6, pi/2, -pi/2, 5*pi/6, -5*pi/6,...
    0, pi/3, -pi/3, 2*pi/3, -2*pi/3, pi,...
    pi/6, -pi/6, pi/2, -pi/2, 5*pi/6, -5*pi/6 ];
% di = sqrt(ri.^2+r^2-2*ri*r.*cos(phi-theta));

V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 5;
mu = 2;
frac = (hb-hr)/(ht-hr);
temp = 2/pi*V*frac;
lamB = 0.1;
C = temp*lamB;

% condProb = @(r,theta,D) r .*cos(theta)

condProb =  @(r,theta,D) prod((C/mu.*sqrt(ri.^2+r.^2-2*ri.*r.*cos(phi-theta)))./...
    (1+C/mu.*sqrt(ri.^2+r.^2-2*ri.*r.*cos(phi-theta)))).*2.*r/D^2 * 1/(2*pi);

int(int(condProb(r,theta,D), 0,2*pi), 0, D) 

