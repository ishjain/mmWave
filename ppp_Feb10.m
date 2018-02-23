%PPP
%Feb9: Verify the theoretical vs emperical model of E(prod(di))

clear all;
clc;
close all;
% lambda=.10;
% %Method 1
% pproc  = poissrnd(lambda, 100, 2);
% size(pproc)
% plot(pproc(:, 1), pproc(:, 2), '.');
% title('Poisson with poissrnd')

rho_b = 0.65;%Rajeev calculated central park
V=1;
hb = 1.8;
hr = 1.4;
ht = 6;
frac = (hb-hr)/(ht-hr);
temp = 2/pi*rho_b*V*frac;
mu = 7;
c = temp/mu;
wannaplot=0;
Rval = 100;%sqrt(10^6/pi);
Lval=(.1:.1:2)/100^2;%(.5:.5:10)/(1000)^2;
count=0;

for iter = 1:100
%     count=count+1;
    for Rind=1:length(Rval)
        for Lind = 1:length(Lval)
            R=Rval(Rind);%m
            
            lambda = Lval(Lind);
            npoints = poissrnd(lambda*pi*R^2);
%             if (npoints==0), count=count+1; continue; end;
            r = R*sqrt(rand(npoints,1));
            
            
            if(wannaplot==1)
                theta = 2*pi*rand(npoints,1);
                x0=0; y0=0; %origin UE
                x = x0 + r.*cos(theta);
                y = y0 + r.*sin(theta);
                scatter(x,y)
                th=0:0.01:2*pi; xx = x0+R*cos(th); yy = y0 + R*sin(th);
                hold on;
                plot(xx,yy)
            end
            %         prod_rsq(iter, Rind) = prod(r.^2);
            %         prod_r(iter, Rind) = prod(r);
            prod_noapprox(iter, Lind) = prod(c*r./(1+c*r));
        end
    end
end
emp_noapprox = mean(prod_noapprox,1);

% emp = mean(prod_rsq,1);
% a_r = 2/3*Rval;
% th_r = exp(-(1-a_r).*pi.*Rval.^2*lambda);
% th_rsq = exp(-(1-Rval.^2/2).*pi.*Rval.^2*lambda);

% a_noapprox = 1-
th_noapprox = exp(-2*pi.*Rval.*Lval/c).*(1+c*Rval).^(2*pi.*Lval/c^2);
% Lval=Lval*
figure(1)
semilogy(Lval*1000^2,emp_noapprox);
hold on;
semilogy(Lval*1000^2, th_noapprox)
legend('emperical', 'theoretical')

h=figure(2);
% plot(Lval,emp_noapprox,'LineWidth',2);
hold on; grid on;
plot(Lval*1000^2, th_noapprox,'LineWidth',2)
% legend('emperical', 'theoretical')
xlabel('\lambda_T (Density of APs per km^2)', 'fontsize',13)
ylabel('Probability that all APs are blocked','fontsize',13)
% title('APs blockage probability in 1km^2 area')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures/N_block_prob.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures/N_block_prob.png'],'-dpng','-r0')



%%Update Feb 20
probNoConnection = exp(-Lval*pi*R^2);%basically n=0 probability
probBl_givenConnection = th_noapprox - probNoConnection;

h=figure(3);
hold on; grid on;
plot(Lval*1000^2, probBl_givenConnection,'LineWidth',2)
xlabel('\lambda_T (Density of APs per km^2)', 'fontsize',13)
ylabel('Prob all-blocked given atleast 1 AP in 100m','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[pwd '/figures/N_block_prob_cond.pdf'],'-dpdf','-r0')
print(h,[pwd '/figures/N_block_prob_cond.png'],'-dpng','-r0')

