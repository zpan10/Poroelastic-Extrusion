% Read me: this program runs fibers_x_nmr.m and plots the figures in
% Figure 2. The displacement plot uses computed data at different time
% steps.

tic

clear all; close all; clc;

dir_prefix = './figures/';

spongepar.plot_switch = 0;

% Euler coordinates of interest
spongepar.x1 = 0.4;
spongepar.x2 = 0.7;

% initial porosity
n0 = 0.85;
spongepar.n0 = n0;

% barrel nozzle area ratio chi
r = 5;  %r(1/n0)<1)
spongepar.r = r;

spongepar.bc_fluidx = 'q-fixed'; % flow rate fixed 

spongepar.dp = NaN;
spongepar.use_force=-1; % ignore existing result

% qall = v_a
qall = 0.63;
spongepar.q = qall; 

% time steps
tall = zeros(length(qall),1);
tobs = [0,linspace(10^(-5),0.5317,200)];
spongepar.tobs = tobs;


spongepar.perm_law = 'Fb';
spongepar.stress_law = 'linear';
[xss,ts,nss,uss,pss,sigss,qs,dps,deltas,v1,v2,sig1] = fibers_x_nmr(dir_prefix,spongepar);
    
      
%% Plot delta
figure
figw = 17.2;%cm
figh = 0.75*figw;
plot(ts,deltas,'-k','LineWidth',4)

ylabel('$\delta/L$','Interpreter','latex')
% ylabel('$\phi_f$','Interpreter','latex')
set(gca,'XLim',[0,1],'YLim',[0,1])
set(gca,'FontSize',24,'FontName','Times New Roman','LineWidth',3)

xlabel('$t/T_{pe}$','Interpreter','latex')

set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 figw figh])
set(gcf,'PaperSize',[figw figh])


%% Plot velocities at blue and red locations
figure
figw = 17.2;%cm
figh = 0.75*figw;
box on
hold on

v1_=v1/qall/r;
v2_=v2/qall/r;

plot(tobs,v1_,'-','LineWidth',4,'Color',[1 .5 .5])
plot(tobs,v2_,'-','LineWidth',4,'Color',[.5 .5 1])

ylabel('$v_s/\chi v_a$','Interpreter','latex')
% ylabel('$\phi_f$','Interpreter','latex')
set(gca,'XLim',[0,0.51],'YLim',[0,1.2])
set(gca,'FontSize',24,'FontName','Times New Roman','LineWidth',3)

xlabel('$t/T_{pe}$','Interpreter','latex')

set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 figw figh])
set(gcf,'PaperSize',[figw figh])


%% Plot displacement
figure
CT=cbrewer('seq', 'Blues', 9);
colormap(CT)
figw = 17.2;%cm
figh = figw;
% t(i) i = 10, 40, 70, 100, 130, 160, 190

t=[10,40,70,100,130,160,190];

delta = ts.*qall;

box on
hold on
plot([delta(10),1],[delta(10),delta(10)],'--','Color',CT(end-6,:),'LineWidth',2)
plot([delta(40),1],[delta(40),delta(40)],'--','Color',CT(end-5,:),'LineWidth',2)
plot([delta(70),1],[delta(70),delta(70)],'--','Color',CT(end-4,:),'LineWidth',2)
plot([delta(100),1],[delta(100),delta(100)],'--','Color',CT(end-3,:),'LineWidth',2)
plot([delta(130),1],[delta(130),delta(130)],'--','Color',CT(end-2,:),'LineWidth',2)
plot([delta(160),1],[delta(160),delta(160)],'--','Color',CT(end-1,:),'LineWidth',2)
plot([delta(190),1],[delta(190),delta(190)],'--','Color',CT(end,:),'LineWidth',2)

plot(xss(10,1:end),uss(10,1:end),'-','Color',CT(end-6,:),'LineWidth',3)
plot(xss(40,1:end),uss(40,1:end),'-','Color',CT(end-5,:),'LineWidth',3)
plot(xss(70,1:end),uss(70,1:end),'-','Color',CT(end-4,:),'LineWidth',3)
plot(xss(100,1:end),uss(100,1:end),'-','Color',CT(end-3,:),'LineWidth',3)
plot(xss(130,1:end),uss(130,1:end),'-','Color',CT(end-2,:),'LineWidth',3)
plot(xss(160,1:end),uss(160,1:end),'-','Color',CT(end-1,:),'LineWidth',3)
plot(xss(190,1:end),uss(190,1:end),'-','Color',CT(end,:),'LineWidth',3)

xlabel('$x/L$','Interpreter','latex')
ylabel('$u_s/L$','Interpreter','latex')
% ylabel('$\phi_f$','Interpreter','latex')
set(gca,'XLim',[0,1],'YLim',[0,1])
set(gca,'FontSize',24,'FontName','Times New Roman','LineWidth',3)
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 figw figh])
set(gcf,'PaperSize',[figw figh])

%% Normalized elastic stress
figure
CT=cbrewer('seq', 'Blues', 9);
colormap(CT)
figw = 17.2;%cm
figh = figw;
% t(i) i = 10, 40, 70, 100, 130, 160, 190

subplot(211)
box on
hold on
plot(xss(10,1:end),sigss(10,1:end),'-','Color',CT(end-6,:),'LineWidth',3)
plot(xss(40,1:end),sigss(40,1:end),'-','Color',CT(end-5,:),'LineWidth',3)
plot(xss(70,1:end),sigss(70,1:end),'-','Color',CT(end-4,:),'LineWidth',3)
ylabel('|\sigma^{\prime}_{xx}/{\it E}_{eff}|')
% ylabel('$\phi_f$','Interpreter','latex')
set(gca,'XLim',[0,1],'YLim',[0,0.6]) % sig
% set(gca,'XLim',[0,1],'YLim',[0.85,1]) % phi
set(gca,'FontSize',24,'FontName','Times New Roman','LineWidth',3)

subplot(212)
box on
hold on
plot(xss(100,1:end),sigss(100,1:end),'-','Color',CT(end-3,:),'LineWidth',3)
plot(xss(130,1:end),sigss(130,1:end),'-','Color',CT(end-2,:),'LineWidth',3)
plot(xss(160,1:end),sigss(160,1:end),'-','Color',CT(end-1,:),'LineWidth',3)
plot(xss(190,1:end),sigss(190,1:end),'-','Color',CT(end,:),'LineWidth',3)
xlabel('$x/L$','Interpreter','latex')
ylabel('|\sigma^{\prime}_{xx}/{\it E}_{eff}|')
% ylabel('$\phi_f$','Interpreter','latex')
set(gca,'XLim',[0,1],'YLim',[0,0.6]) % sig
% set(gca,'XLim',[0,1],'YLim',[0.85,1]) %phi
set(gca,'FontSize',24,'FontName','Times New Roman','LineWidth',3)
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 figw figh])
set(gcf,'PaperSize',[figw figh])




