% Read me: this program runs fibers_x_nmr.m at various v_a from 0.01 to 10
% and return the extrusion time tall (scaled by Tpe). The user needs to 
% define barrel nozzle area ratio r and initial porosity n0 for each run.
% 
% Usual run time ~1h
% 
% User set value:
% For figure 3a 
% n0: 0.95, r: 5, 20, 10 
% For figure 3d 
% n0:0.9 r: 10,8,6,4,2l

tic

clear all; close all; clc;

dir_prefix = './figures/';

spongepar.plot_switch = 0;

% initial porosity
n0 = 0.95;
spongepar.n0 = n0;

% barrel nozzle area ratio chi
r = 5;  %r(1/n0)<1)
% r = 20;
% r = 10;
spongepar.r = r;

spongepar.bc_fluidx = 'q-fixed'; % flow rate fixed 

spongepar.dp = NaN;
spongepar.use_force=-1; % ignore existing result

% dimensionless v_a 
qall = (logspace(log10(0.01),log10(10),10))';
spongepar.q = qall; 

% time steps
tall = zeros(length(qall),1);
tobs = [0,linspace(10^(-5),0.02,200)];
spongepar.tobs = tobs;


% Loop over different v_a
for l=1:length(qall)
    spongepar.q = qall(l); 
i0 = 1/qall(l)/r*0.9;
i = i0;

% Find the finishing time when deltas(end) > 0.99
    while 1
    i = i+0.001;
    tobs = [0,linspace(10^(-5),i,200)];
    spongepar.tobs = tobs;
    
    spongepar.perm_law = 'Fb';
    spongepar.stress_law = 'linear';
    [xss,ts,nss,uss,pss,sigss,qs,dps,deltas,v1,sig1] = fibers_x_nmr(dir_prefix,spongepar);
    
        if deltas(end) < 0.7
            i = i + 0.2*i0;
        end

        if deltas(end) < 0.85
            i = i + 0.08*i0;
        end 

        if deltas(end) < 0.98
            i = i + 0.001*i0;
        end 

        tall(l) = i;

        if deltas(end) > 0.99
            toc
            break
        end

    end
end


% ratio for volume fraction increase
phi_ex_0 = 1./qall./tall;

x=logspace(-1,0.1)';
vq = interp1(qall,phi_ex_0,x,'spline');
Chi5=[x vq];


CT=cbrewer('seq', 'Oranges', 20);
colormap(CT)


% semilogx(Chi20(:,1),Chi20(:,2),'-','Color', CT(20,:),'LineWidth',10);
% hold on
% semilogx(Chi12(:,1),Chi12(:,2),'-','Color', CT(12,:),'LineWidth',10);
semilogx(Chi5(:,1),Chi5(:,2),'-','Color', CT(5,:),'LineWidth',10);

figw = 40;%cm
figh = figw*3/4;
xlabel('$\overline{v}_a$','Interpreter','latex','Fontweight','bold')
ylabel('$\frac{\phi_{s,ex}}{\phi_{s,0}}$','Interpreter','latex','Fontweight','bold')
set(gca,'XLim',[0.01,2],'YLim',[0,22])
set(gca,'FontSize',20,'FontName','Times New Roman','LineWidth',4.5)
set(gcf,'PaperUnits','centimeters','Color','none')
set(gcf,'PaperPosition',[0 0 figw figh])
set(gcf,'PaperSize',[figw figh])


