function data = example_cie_target_tracking
% --- example_cie_target_tracking() ---------------------------------------
% Example of the common information estimate (CIE) defined in the paper: 
% "Decentralized Data Fusion of Dimension-Reduced Estimates Using Local 
% Information Only", IEEE Aerospace Conference 2023.  
%
% The example is a decentralized data fusion scenario where two agents, 
% each equipped with a local sensor and a Kalman filter, in 2D measure and 
% track a common target. The agents communicate dimension-reduced estimate 
% derived using GEVO-KF and GEVO-CI as described in the abovementioned 
% paper. KF and CI runs in parallel and uses the same measurements. 
%
% The purpose with this example is to further illustrate the CIE using code
% and visualizations. The simulation is visualized using plots of the local 
% estimates and plots of the J(P) and J(Gamma), where P is the local 
% estimate covariance, Gamma (G) is the common information estimate 
% covariance and J() = trace(). Especially the plots with J(Gamma) give
% insights into how the CIE is updated with respect to the local estimate.
%
% The communication is according to the basic communication schedule.
% 
% 2023-10-30 Robin Forsling

set_latex_interpreter

close all


% --- PARAMETERS ----------------------------------------------------------
par = get_parameters;
par.m = 2;


% --- TARGET --------------------------------------------------------------
XT = get_circle_arc_fixed_length('right',[0;0],-pi/6,par.r,par.L,par.N);


% --- AGENTS --------------------------------------------------------------
[RS1,RS2] = get_measurement_covariances;
agent = cell(2,1);
agent{1} = get_agent(par,RS1);
agent{2} = get_agent(par,RS2);


% --- PLOT SETTINGS -------------------------------------------------------
plot_par = get_plot_parameters(par,XT,agent);


% --- RUN SIMULATIONS -----------------------------------------------------
C = cell(2,1); data.J = C;
data.xhat_kf = C; data.P_kf = C; 
data.xhat_ci = C; data.P_ci = C;
for i = 1:2
    A = zeros(3,par.N);                                                     % row 1: before communication, row 2: after communication before fusion, row 3 : after fusion
    data.J{i}.P_kf = A; data.J{i}.G_kf = A;
    data.J{i}.P_ci = A; data.J{i}.G_ci = A;
    A = zeros(par.nx,par.N); B = zeros(par.nx,par.nx,par.N);
    data.xhat_kf{i} = A; data.P_kf{i} = B; 
    data.xhat_ci{i} = A; data.P_ci{i} = B; 
end

y = cell(2,1);
for k = 1:par.N
    

    % --- 0. GENERATE MEAUREMENTS -----------------------------------------
    for i = 1:2; y{i} = generate_measurement(XT(:,k),agent{i}); end


    % --- 1. INITIALIZE USING FIRST MEASUREMENT ---------------------------
    if k == 1
        for i = 1:2
            agent{i}.kf.xhat = zero_pad_end(y{i},agent{i}.pm.nx);
            agent{i}.kf.P = agent{i}.pm.P0;
            agent{i}.kf.ghat = agent{i}.kf.xhat;
            agent{i}.kf.G = agent{i}.pm.G0;
            agent{i}.ci = agent{i}.kf;
        end


    % --- 2. PREDICT AND UDPATE LOCAL ESTIMATES ---------------------------
    else
        for i = 1:2
            agent{i} = predict_estimates(agent{i});
            agent{i} = update_estimates(agent{i},y{i});
        end
    end

    data = save_data(1,data,agent,k);


    % --- 3. RUN GEVO AND SIMULATE COMMUNICATION --------------------------
    [itx,irx] = get_communication_topology(k);
    [agent{itx},tx_est] = run_gevo(par,agent{itx});

    data = save_data(2,data,agent,k);


    % --- 4. FUSE RECEIVED ESTIMATE ---------------------------------------
    agent{irx} = fuse_received_estimate(agent{irx},tx_est);

    data = save_data(3,data,agent,k);


    % --- PLOT ------------------------------------------------------------
    plot_estimates(plot_par,k,XT,data); 
    plot_loss_functions(plot_par,k,data);

end

end


% --- SENSOR AND FILTER FUNCTIONALITY -------------------------------------
function y = generate_measurement(xt,agent)
    L = chol(agent.sm.R,'lower');
    v = L*randn(agent.sm.ny,1);
    y = xt + v;                                                             % xt only contains position
end

function agent = predict_estimates(agent)
    F = agent.pm.F; Q = agent.pm.Q;

    % KF local estimate:
    agent.kf.xhat = F*agent.kf.xhat;
    agent.kf.P = F*agent.kf.P*F' + Q;

    % KF common information estimate:
    agent.kf.ghat = F*agent.kf.ghat;
    agent.kf.G = F*agent.kf.G*F' + Q;

    % CI local estimate:
    agent.ci.xhat = F*agent.ci.xhat;
    agent.ci.P = F*agent.ci.P*F' + Q;

    % CI common information estimate:
    agent.ci.ghat = F*agent.ci.ghat;
    agent.ci.G = F*agent.ci.G*F' + Q;    
end

function agent = update_estimates(agent,y)
    H = agent.sm.H; R = agent.sm.R; I = eye(agent.pm.nx);

    % KF local estimate:
    ybar = y - H*agent.kf.xhat;
    S = H*agent.kf.P*H' + R;
    K = agent.kf.P*H'/S;
    agent.kf.xhat = agent.kf.xhat + K*ybar;
    agent.kf.P = (I - K*H)*agent.kf.P;

    % CI local estimate:
    ybar = y - H*agent.ci.xhat;
    S = H*agent.ci.P*H' + R;
    K = agent.ci.P*H'/S;
    agent.ci.xhat = agent.ci.xhat + K*ybar;
    agent.ci.P = (I - K*H)*agent.ci.P;
end

function [itx,irx] = get_communication_topology(k)
    if mod(k,2) == 0; itx = 2; irx = 1;
    else; itx = 1; irx = 2;
    end
end

function agent = fuse_received_estimate(agent,rx_est)

    % KF
    [agent.kf.xhat,agent.kf.P] = fuse_kf(agent.kf.xhat,agent.kf.P,rx_est.kf);
    [agent.kf.ghat,agent.kf.G] = fuse_kf(agent.kf.ghat,agent.kf.G,rx_est.kf);

    % CI
    [agent.ci.xhat,agent.ci.P] = fuse_ci(agent.ci.xhat,agent.ci.P,rx_est.ci);
    [agent.ci.ghat,agent.ci.G] = fuse_ci(agent.ci.ghat,agent.ci.G,rx_est.ci);    
end

function [xhat,P] = fuse_kf(y1,R1,est)
    Psi = est.Psi; yPsi = est.yPsi; RPsi = est.RPsi;
    I1 = inv(R1); IPsi = inv(RPsi);
    P = inv(I1 + Psi'*IPsi*Psi);
    xhat = P*(I1*y1 + Psi'*IPsi*yPsi);
end

function [xhat,P] = fuse_ci(y1,R1,est)
    Psi = est.Psi; yPsi = est.yPsi; RPsi = est.RPsi;
    I1 = inv(R1); IPsi = inv(RPsi);
    f = @(w) trace(inv(w*I1 + (1-w)*Psi'*IPsi*Psi));
    w = fminbnd(f,0,1);
    P = inv(w*I1 + (1-w)*Psi'*IPsi*Psi);
    xhat = P*(w*I1*y1 + (1-w)*Psi'*IPsi*yPsi);
end


% --- GEVO FUNCTIONALITY --------------------------------------------------
function [agent,tx_est] = run_gevo(par,agent)
    
    % KF
    y1 = agent.kf.ghat; R1 = agent.kf.G;
    [y2,R2] = subtract_estimate(agent.kf.xhat,agent.kf.P,y1,R1);
    Psi = gevo_kf(R1,R2,par.m);
    tx_est.kf.Psi = Psi;
    tx_est.kf.yPsi = Psi*y2;
    tx_est.kf.RPsi = Psi*R2*Psi';
    [agent.kf.ghat,agent.kf.G] = fuse_kf(y1,R1,tx_est.kf);

    % CI
    y1 = agent.ci.ghat; R1 = agent.ci.G;
    y2 = agent.ci.xhat; R2 = agent.ci.P;
    Psi = gevo_ci(R1,R2,par.m);
    tx_est.ci.Psi = Psi;
    tx_est.ci.yPsi = Psi*y2;
    tx_est.ci.RPsi = Psi*R2*Psi';
    [agent.ci.ghat,agent.ci.G] = fuse_ci(y1,R1,tx_est.ci);
end

function Psi = gevo_kf(R1,R2,m)
    Q = R1*R1;
    S = R1+R2; 
    [X,D] = eig(Q,S,'qz');
    idx_vec = get_max_idx_vec(D,m);
    Psi = gram_schmidt_process(X(:,idx_vec))';
end

function Psi = gevo_ci(R1,R2,m)
    max_iter = 50; J_thres = 0.0001; wmin = 0.01; wmax = 1;                     % Algorithm settings
    I1 = inv(R1);   
    w = 0.5; J = Inf;
    for k = 1:max_iter
        J_prev = J;
        Q = R1*R1/(w*w);
        S = R1/w + R2/(1-w);
    
        % SOLVE FOR Psi
        [X,D] = eig(Q,S,'qz');
        idx_vec = get_max_idx_vec(D,m);
        Psi = X(:,idx_vec)';
    
        % SOLVE FOR w
        f = @(w) trace(inv(w*I1 + (1-w)*Psi'/(Psi*R2*Psi')*Psi));
        w = fminbnd(f,wmin,wmax);
        J = ci_loss_function(R1,R2,Psi,w);
    
        if (J_prev-J)/J < J_thres; break; end
    end    
    Psi = gram_schmidt_process(Psi);
end

function J = ci_loss_function(R1,R2,Psi,w)
    if w >= 1; J = trace(R1); return;
    elseif w <= 0; error('ci_loss_function: omega out of range')
    end
    J = trace(inv(w*inv(R1) + (1-w)*Psi'/(Psi*R2*Psi')*Psi));
end

function U = gram_schmidt_process(X)
    [n,m] = size(X);
    is_tall = n > m;
    if ~is_tall; X = X'; [n,m] = size(X); end                         
    U = zeros(n,m);
    U(:,1) = X(:,1)/norm(X(:,1));
    for i = 2:m
        x = X(:,i); u = x;
        for j = 1:i-1
            uprev = U(:,j);
            u = u - projection_operator(uprev,x);
        end
        U(:,i) = u/norm(u);
    end
    if ~is_tall; U = U'; end 
end

function vu = projection_operator(u,v)
    vu = (u'*v/(u'*u))*u;
end

function [c,C] = subtract_estimate(a,A,b,B)
    Ai = inv(A); Bi = inv(B);
    C = inv(Ai - Bi);
    c = C*(Ai*a - Bi*b);
end


% --- MOTION --------------------------------------------------------------
function x = get_circle_arc_fixed_length(turn_dir,x0,head0,r,L,N)
    x = zeros(2,N); 
    turn_ang = L/r;
    ang_vec = turn_ang*linspace(0,N,N)/N;
    for k = 2:N
        ang = ang_vec(k);
        if strcmpi(turn_dir,'left'); x(:,k) = r*[sin(ang) ; (1-cos(ang))];  % Left hand turn
        else; x(:,k) = r*[sin(ang) ; (cos(ang)-1)];                         % Right hand turn
        end
    end
    x = rotate_motion(x,[0;0],head0) + x0(1:2);
end

function xrot = rotate_motion(x,xc,alpha)
    T = [cos(alpha) -sin(alpha) ; sin(alpha) cos(alpha)];
    xrot = T*(x-xc)+xc;
end


% --- PLOT FUNCTIONALITY --------------------------------------------------
function plot_estimates(plot_par,k,XT,data) 
    f = figure(1);clf;
    f.Position = [100 100 400 1200];
    for i = 1:2
        xkf = data.xhat_kf{i}; xci = data.xhat_ci{i};
        subplot(2,1,i);hold on
        h = plot(XT(1,1:k),XT(2,1:k),'k--','DisplayName','Target trajectory'); h.LineWidth = 1.5; 
        h = plot(XT(1,k),XT(2,k),'ko','DisplayName','Target position'); h.MarkerFaceColor = [0 0 0]; h.MarkerSize = 5;
        h = plot(xkf(1,1:k),xkf(2,1:k),'-','DisplayName','KF estimated trajectory'); h.LineWidth = plot_par.lw; h.Color = plot_par.c.kf;
        h = plot(xkf(1,k),xkf(2,k),'*','DisplayName','KF estimated position'); h.MarkerSize = plot_par.ms; h.Color = plot_par.c.kf;
        h = plot(xci(1,1:k),xci(2,1:k),'-','DisplayName','CI estimated trajectory'); h.LineWidth = plot_par.lw; h.Color = plot_par.c.ci;
        h = plot(xci(1,k),xci(2,k),'x','DisplayName','CI estimated position'); h.MarkerSize = plot_par.ms; h.Color = plot_par.c.ci;
        xlabel('$x_1$','interpreter','latex'); ylabel('$x_2$','interpreter','latex'); title(sprintf('Agent %d: Estimate',i))
        box on; %axis equal
        xlim(plot_par.xlim2d); ylim(plot_par.ylim2d)
        legend('show','location','southwest')  

        set_fontsize_all(14)
    end
end

function plot_loss_functions(plot_par,k,data)
    t = 1:k;
    f = figure(2);clf;
    f.Position = [600 100 800 1200];
    for i = 1:2
        subplot(2,2,i);hold on
        P = data.J{i}.P_kf; G = data.J{i}.G_kf;
        h = plot(t,P(1,1:k),':','DisplayName','$J(P)$ Before communication'); h.Color = plot_par.c.P; h.LineWidth = plot_par.lw;
        h = plot(t,P(2,1:k),'--','DisplayName','$J(P)$ After comm., before fus.'); h.Color = plot_par.c.P; h.LineWidth = plot_par.lw;
        h = plot(t,P(3,1:k),'-','DisplayName','$J(P)$ After fusion'); h.Color = plot_par.c.P; h.LineWidth = plot_par.lw;
        h = plot(t,G(1,1:k),':','DisplayName','$J(\Gamma)$ Before communication'); h.Color = plot_par.c.G; h.LineWidth = plot_par.lw;
        h = plot(t,G(2,1:k),'--','DisplayName','$J(\Gamma)$ After comm., before fus.'); h.Color = plot_par.c.G; h.LineWidth = plot_par.lw;
        h = plot(t,G(3,1:k),'-','DisplayName','$J(\Gamma)$ After fusion'); h.Color = plot_par.c.G; h.LineWidth = plot_par.lw;        

        xlabel('$k$','interpreter','latex'); ylabel('$J$','interpreter','latex'); title(sprintf('Agent %d: KF',i))
        box on; grid on
        xlim(plot_par.xlimJ); ylim(plot_par.ylimJ)
        legend('show','NumColumns',1,'location','northeast')

        subplot(2,2,i+2);hold on
        P = data.J{i}.P_ci; G = data.J{i}.G_ci;
        h = plot(t,P(1,1:k),':','DisplayName','$J(P)$ Before communication'); h.Color = plot_par.c.P; h.LineWidth = plot_par.lw;
        h = plot(t,P(2,1:k),'--','DisplayName','$J(P)$ After comm., before fus.'); h.Color = plot_par.c.P; h.LineWidth = plot_par.lw;
        h = plot(t,P(3,1:k),'-','DisplayName','$J(P)$ After fus.'); h.Color = plot_par.c.P; h.LineWidth = plot_par.lw;
        h = plot(t,G(1,1:k),':','DisplayName','$J(\Gamma)$ Before communication'); h.Color = plot_par.c.G; h.LineWidth = plot_par.lw;
        h = plot(t,G(2,1:k),'--','DisplayName','$J(\Gamma)$ After comm., before fus.'); h.Color = plot_par.c.G; h.LineWidth = plot_par.lw;
        h = plot(t,G(3,1:k),'-','DisplayName','$J(\Gamma)$ After fusion'); h.Color = plot_par.c.G; h.LineWidth = plot_par.lw;        

        xlabel('$k$','interpreter','latex'); ylabel('$J$','interpreter','latex'); title(sprintf('Agent %d: CI',i))
        box on; grid on
        xlim(plot_par.xlimJ); ylim(plot_par.ylimJ)
        %legend('show','NumColumns',1,'location','northeast')
        set_fontsize_all(14)
    end
end

function plot_par = get_plot_parameters(par,XT,agent)

    xmin = min(XT(1,:)); xmax = max(XT(1,:)); xc = 0.5*(xmax+xmin);
    ymin = min(XT(2,:)); ymax = max(XT(2,:)); yc = 0.5*(ymax+ymin); 
    L = 1.2*max(abs([(xmax-xmin) (ymax-ymin)]));
    plot_par.xlim2d = [-L L] + xc;
    plot_par.ylim2d = [-L L] + yc;
    plot_par.xlimJ = [0 par.N+1];
    plot_par.ylimJ = [0 1.2*trace(agent{1}.pm.G0)];

    plot_par.lw = 1.0;
    plot_par.ms = 7;
    plot_par.c = get_colors;
end

function c = get_colors()
    clr = get_thesis_colors;
    c.kf = clr.kf;
    c.ci = clr.ci;
    c.P = clr.darkyellow;
    c.G = clr.blue;
end


% --- MISC FUNCTIONALITY --------------------------------------------------
function data = save_data(irow,data,agent,k)
    for i = 1:2 % SAVE DATA
        data.J{i}.P_kf(irow,k) = trace(agent{i}.kf.P); 
        data.J{i}.G_kf(irow,k) = trace(agent{i}.kf.G); 
        data.J{i}.P_ci(irow,k) = trace(agent{i}.ci.P); 
        data.J{i}.G_ci(irow,k) = trace(agent{i}.ci.G); 
        if irow == 3
            data.xhat_kf{i}(:,k) = agent{i}.kf.xhat; 
            data.P_kf{i}(:,:,k) = agent{i}.kf.P;
            data.xhat_ci{i}(:,k) = agent{i}.ci.xhat; 
            data.P_ci{i}(:,:,k) = agent{i}.ci.P;
        end
    end
end

function par = get_parameters()
    par.m = 1;
    par.ny = 2;                                                             % Measurement dimensionality, should be 2
    par.nx = 6;                                                             % State dimensionality, should be 6
    par.N = 12;                                                             % Number of time steps
    par.r = 200;
    par.L = 120;
    par.q = 1;                                                              % Process noise magnitude
    par.Ts = 1;                                                             % Sampling time
    par.vmax = 10;
    par.amax = 1;
end

function a = get_agent(par,RS)   
    Ts = par.Ts; ny = par.ny; nx = par.nx; 

    a.pm.nx = nx;
    a.pm.q = par.q;
    a.pm.Ts = Ts;
    a.pm.F = [  eye(ny) Ts*eye(ny) 0.5*Ts^2*eye(2) ; ...
                zeros(ny) eye(ny) Ts*eye(ny) ; ...
                zeros(ny) zeros(ny) eye(ny)];
    a.pm.Q = par.q^2* ...
             [  Ts^5*eye(ny)/20 Ts^4*eye(ny)/8 Ts^3*eye(ny)/6 ; ...
                Ts^4*eye(ny)/8 Ts^3*eye(ny)/3 Ts^2*eye(ny)/2 ; ...
                Ts^3*eye(ny)/6 Ts^2*eye(ny)/2 Ts*eye(ny)];
    a.pm.P0 = blkdiag(RS,par.vmax^2*eye(ny),par.amax^2*eye(ny));
    a.pm.G0 = 2*a.pm.P0;

    a.sm.ny = ny;
    a.sm.H = eye(ny,nx);                       
    a.sm.R = RS;

    a.kf.xhat = [];
    a.kf.P = [];
    a.kf.ghat = [];
    a.kf.G = [];

    a.ci = a.kf;
end

function idx_vec = get_max_idx_vec(val,m)   
    if ~isvector(val); val = diag(val); end
    idx_vec = zeros(1,m);
    for k = 1:m
        [~,idx] = max(val);
        idx_vec(k) = idx;
        val(idx) = -Inf;
    end
end

function [R1,R2] = get_measurement_covariances()
    T = @(a) [cos(a) -sin(a) ; sin(a) cos(a)];
    a1 = pi/9; a2 = -pi/9;
    R1 = T(a1)*diag([16 1])*T(a1)';
    R2 = T(a2)*diag([16 1])*T(a2)';
end

function v = zero_pad_end(v,n)
    d = n-length(v);
    v = [v ; zeros(d,1)];
end


