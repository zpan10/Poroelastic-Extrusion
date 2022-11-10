function [xss,ts,nss,uss,pss,sigss,qs,dps,as,vx1,vx2,sigx1] = fibers_x_nmr(dir_prefix,spongepar)
    % Scaling:
    %   n is porosity phi_f
    %   x is x/L
    %   t is t/T
    %   u is u/L: displacement 
    %   p is p/M
    %   r is vnoz/vavg
    %   sig is sig/M: stress
    %   q is mu*q*L/(k0*M)
    %   dp is dp/M = (p(delta)-p(L))/M
    %   v is v/L*T = mu*v*L/(k0*M)
    %   a is a/L
    %   adot is (T/L)*adot = (T/L)*da/dt
    %
    % Where
    %   n0 is the initial porosity
    %   L is the initial length of the sponge
    %   T = mu*(L^2)/(k0*M);
    %
    % Geometry:
    %   This code takes the left and right edges of the sponge to be at x=a(t) and x=b=L, respectively,
    %   where the left edge is free and the right edge is fixed. As in the paper, the code assumes 
    %   without loss of generality that a0/L=0 and b0/L=b/L=1.
    %
    % INPUTS
    %   dir_prefix: Directory in which to save the results.
    %   spongepar: Structure with parameters and options.
    %     *.q: Dimensionless flow rate
    %     *.dp: Dimensionless pressure drop
    %       - Exactly one of the two above must be passed in. The other should be absent or NaN.
    %     *.sigstar: Dimensionless effective stress applied at x=a(t)
    %     *.n0: Initial (relaxed) porosity
    %     *.tobs: Times at which to save and return the solution.
    %
    % OPTIONAL INPUTS
    %   spongepar
    %     *.plot_switch: Plot occasionally while solving [1, default], or not [0].
    %     *.use_force: If there is an existing results file, ignore it [-1] or use it [any other value, default]
    %     *.stress_law: 
    %       - linear: Linear elasticity [default]
    %       - log: Hencky elasticity
    %     *.perm_law:
    %       - const: Constant permeability, where k is k/k0 = 1 [default]
    %       - KC: normalized Kozeny-Carman, where k is k/k0 = (((1-n0)^2)/(n0^3))*((n^3)/((1-n)^2))
    
    % Set some default values
    if ~isfield(spongepar,'use_force')
       spongepar.use_force = 0;
    end
    use_force = spongepar.use_force;
    if ~isfield(spongepar,'plot_switch')
        spongepar.plot_switch = 0;
    end
    plot_switch = spongepar.plot_switch;
    if ~isfield(spongepar,'stress_law')
        spongepar.stress_law = 'linear';
    end
    stress_law = spongepar.stress_law;
    if ~isfield(spongepar,'perm_law')
        spongepar.perm_law = 'const';
    end
    perm_law = spongepar.perm_law;
    
    % Unpack physical parameters
    n0 = spongepar.n0; % [-] initial porosity
    x1 = spongepar.x1; % Euler location of interest
    x2 = spongepar.x2; % Euler location of interest
    if ~isfield(spongepar,'q')  % specify q -ZP
        spongepar.q = NaN;
    end
    q = spongepar.q;
    if ~isfield(spongepar,'dp')
        spongepar.dp = NaN;
    end
    dp = spongepar.dp;
    if isnan(q) && ~isnan(dp)
        bc_fluid = 'dp-fixed';
    elseif ~isnan(q) && isnan(dp)
        bc_fluid = 'q-fixed';
    else
        error('Must input exactly ONE of q and dp. The other will be determined by the solution.')
    end
    if ~isfield(spongepar,'sigstar')
        spongepar.sigstar = 0;
    end
    sigstar = spongepar.sigstar;
    
    use = 'n';
    save_filename = ['sponge_x_nmr' '_sig_' stress_law '_k_' perm_law '_q' num2str(q) '_dp' num2str(dp) '_sigstar' num2str(sigstar) '_n0' num2str(n0) '.mat'];
    if exist(strcat(dir_prefix,save_filename))==2 && use_force~=-1
        % Use saved results automatically, unless use_force==-1, in which case ignore them.
        disp('Using existing sponge_nmr results...')
        load(strcat(dir_prefix,save_filename),'spongepar','xss','ts','nss','uss','pss','sigss','qs','dps','as');
    else
        
        % Unpack control parameters
        tobs = spongepar.tobs;
        r = spongepar.r;
        a0 = 0;
		b0 = 1;
		
        N = 400;
        
		% Initial condition: Start in a relaxed state
        a = a0;
        b = b0;
        dx = (b-a)/N;
        xs = [(a+dx/2):dx:(b-dx/2)];
		n0s = n0*ones(1,N);
        
        ns = n0s;
		Vs0 = sum(dx.*(1-n0s)); % For n = porosity
        
        % Define a function for the permeability law
        if strcmp(perm_law,'const')
            k = @(n) ones(size(n)); % Constant
        elseif strcmp(perm_law,'Fb')
            k = @(n) ((1-n0)/(-log(1-n0)-0.931) * ((-log(1-n)-0.931)./(1-n))); % Normalized Fiber 
        else
            error('Unknown permeability law.')
        end
        
        % Define a function sig(n) for the stress law, and its inverse n_from_sig(sig).
        if strcmp(stress_law,'linear')
            sig = @(n) (n-n0)/(1-n0); % Linear elasticity
            n_from_sig = @(s) n0 + (1-n0)*s;
        elseif strcmp(stress_law,'log')
            sig = @(n) log((1-n0)./(1-n))./((1-n0)./(1-n)); % Hencky elasticity
%             n_from_sig = @(s) 1 - (1-n0)*((-s)./lambertw(-s)); % y=-lambertw(-x)/x solves x=ln(y)/y
        else
            error('Unknown stress law.')
        end
        
        % Inner boundary condition: Imposed effective stress
        if sigstar==0
            nstar = n0; % Relaxed. Below should give same result, but this is exact.
        else
            nstar = n_from_sig(sigstar);
        end
        kstar = k(nstar);
        
        if ~isnan(dp)
            % If the pressure is imposed, n(x=b) is fixed in time
            n_at_b = n_from_sig(sigstar-dp);
        else
            % If the flux is imposed n(x=b) changes in time
            n_at_b = NaN;
        end
        
        % Empty structures for saving
        ts = [];
        as = [];
        qs = [];
        dps = [];
        as_anl = [];
        qs_anl = [];
        dps_anl = [];
        
        if plot_switch==1
            fig1 = figure();
        end
        count = 0;
        
        % Initial condition
        Y0 = [(ns');a]; %[0.85; 0.85;...;0.85;0]
        
        options = odeset('RelTol',1E-10, ...
                           'AbsTol',1E-10); 
        %[T,Y,te,ye,ie] = ode45(@odefun,tobs,Y0,options);        
        
        %ode45            
%         [T,Y] = ode45(@odefun,tobs,Y0,options); 
        
        %ode15s
        [T,Y] = ode15s(@odefun,tobs,Y0,options);
        
        
        ts = T;
        nss = Y(:,1:N); % porosity at all times
        as = Y(:,N+1); % left boundary location
        
        [xss,uss,sigss,pss,kss,qs,dps,vx1,vx2,sigx1] = ode_post(nss,as,x1);
        
        % Save the results
        if exist(dir_prefix)~=7
            mkdir(dir_prefix)
        end
        save(strcat(dir_prefix,save_filename),'spongepar','xss','ts','nss','uss','pss','sigss','qs','dps','as','vx1');
        
        if plot_switch==1 && exist('fig1')==1
            close(fig1)
            drawnow()
        end
        
    end
    
    % -----------------------------------------------------
    % -----------------------------------------------------
    
    function [ndots,adot] = calcFF(ns,a)
        
        % Calculate grid
        dx = (b-a)/N;
        xs = (a+dx/2):dx:(b-dx/2);
        
        % Calculate stress field from porosity field
        sigs = sig(ns);
        
        % Calculate permeability field from porosity field
        ks = k(ns);
        kstar = k(nstar); % nstar == n0
        
        % q = q;
        
        adot = q + (2*ks(1)*kstar/(ks(1)+kstar))*(sigs(1)-sigstar)/(dx/2); % sigstar = 0 (sigstar is imposed stress) -ZP
        % adot = q + kstar*(sigs(1)-sigstar)/(dx/2); 
        
        % Transmissibilities (transmissivities?)
        TTs = 2* (1-ns(1:end-1)).*ks(1:end-1).* (1-ns(2:end)).*ks(2:end) ...
            ./(  (1-ns(1:end-1)).*ks(1:end-1) + (1-ns(2:end)).*ks(2:end) );
%         TTs = ( (1-ns(1:end-1)).*ks(1:end-1) + (1-ns(2:end)).*ks(2:end) )/2;
        
        % Fluxes: phi_f*v_f
        Fs_adv = zeros(1,N+1);
        Fs_adv(2:N) = ( q - (((b-(xs(1:end-1)+dx/2))*adot)/(b-a)) ).*((ns(1:end-1)+ns(2:end))/2); % Advective, but not upwinded. More accurate??
%         Fs_adv(1,2:N) = ( q - (((b-(xs(1:end-1)+dx/2))*adot)/(b-a)) ).*ns(1:end-1); 
        Fs_dff = zeros(1,N+1);
        Fs_dff(2:N) = -TTs.*( (sigs(2:end)-sigs(1:end-1))/dx ); % Diffusive: Centered
        Fs = Fs_adv + Fs_dff; % Total
        FL = q-adot;
        FR = q-(1-ns(end))*r*q;  %-ZP
        Fs(1) = FL; Fs(N+1) = Fs(N)+(FR-Fs(N)); % Boundary conditions
        
        
        ndots = (adot/(b-a))*ns - (Fs(1,2:end)-Fs(1,1:end-1))/dx;
        
%         fprintf('%d\n',ndots(N));
%         if ndots(N) <= 1.364880e-01
%           keyboard
%         end
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    function [q] = calc_q(q,sigs,ks)
        
        % !! Leaves q untouched if bc_fluid is q-fixed.
        
        if strcmp(bc_fluid,'dp-fixed')
            
            % Calculate q (needed to evolve n)
            sig_b = sigstar - dp; % p(b)-p(a)=sig(b)-sig(a), dp=p(a)-p(b), sig(a)=sigstar --> sig_b = sigstar-dp
            n_b = n_from_sig(sig_b);
            k_b = k(n_b);
            q = -(2*ks(end)*k_b/(ks(end)+k_b))*(sig_b-sigs(end))/(dx/2);
            % q = -k_b*(sig_b-sigs(end))/(dx/2);
            
        end
        
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    function [dp] = calc_dp(dp,us,sigs)
        
        % !! Leaves dp untouched if bc_fluid is dp-fixed.
        
        if strcmp(bc_fluid,'q-fixed')
            
            % Calculate dp (not needed -- just for reference)
            n_b = n0 + (1-n0)*(us(end)-us(end-1))/(dx/2); % (n-n0)/(1-n0) = du/dx
            sig_b = sig(n_b); 
            dp = sigstar - sig_b; % p(b)-p(a)=sig(b)-sig(a), dp=p(a)-p(b), sig(a)=sigstar --> sig_b = sigstar-dp
            
        end
        
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    function [a_from_u,us] = u_from_n(xs,ns)
        % Calculate the displacement field from the porosity
        
        % (J-1)/J = (n-n0)/(1-n0) = du/dx, where J = (1-n0)/(1-n) is the Jacobian determinant
        % Integrate to get u at the i-1/2 interfaces
        us =  a + [0,cumsum(dx*((ns(1:end)-n0)/(1-n0)))];
        % Add the left boundary value 
        %us = [a,us]; % -ZP
        
        % Save the left boundary value
        a_from_u = a0 + us(1);
        
        % Interpolate to get the cell-centered displacements
        us = interp1([xs-dx/2,xs(end)+dx/2], us, xs);
        
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    function Ydot = odefun(T,Y)
        
        ns = Y(1:N)';
        a = Y(N+1);
        [ndots,adot] = calcFF(ns,a);
        Ydot = [(ndots');adot];
        
    end
    
    function [xss,uss,sigss,pss,kss,qs,dps,vx1,vx2,sigx1] = ode_post(nss,as,x1)
        
        % Empty structures for saving
        xss = zeros(size(nss));
        uss = zeros(size(nss));
        sigss = zeros(size(nss));
        pss = zeros(size(nss));
        kss = zeros(size(nss));
        qs = zeros(size(as));
        dps = zeros(size(as));
        vss = zeros(size(nss));
        vx1 = zeros(size(nss,1),1);
        vx2 = zeros(size(nss,1),1);
        sigx1 = zeros(size(nss,1),1);
        
        % Loop over time
        for it=1:size(nss,1)
            
            ns = nss(it,:);
            a = as(it);
            
            % Calculate grid
            dx = (b-a)/N;
            xs = [(a+dx/2):dx:(b-dx/2)];
            
            % Calculate displacement field from porosity field
            [a_from_u,us] = u_from_n(xs,ns);
            
            % Calculate stress field from porosity field
            sigs = sig(ns);
            sigx1(it) = interp1(xs,sigs,x1);
            
            % Calculate permeability field from porosity field
            ks = k(ns);
            
            v = q + ks(2:end).*(sigs(2:end)-sigs(1:end-1))/dx;
            vss(it,:) = [q + ks(1).*(sigs(2)-sigs(1))/dx, v]; 
            vx1(it) = interp1(xs,vss(it,:),x1);
            
            vx2(it) = interp1(xs,vss(it,:),x2);
            
%             q = calc_q(q,sigs,ks);
%             dp = calc_dp(dp,us,sigs);
            
            % Save
            xss(it,:) = xs;
            uss(it,:) = us;
            sigss(it,:) = sigs;
            %pss(it,:) = (sigs-sigstar) + dp; % dp/dx = dsig/dx, dp = p(a)-p(b)
            pss(it,:) = (sigs-sigstar); % ZP
            kss(it,:) = ks;
            qs(it) = q;
            dps(it) = dp;
            
        end
            
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    
end