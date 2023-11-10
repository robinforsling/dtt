function data = example_cie_convergence
% --- example_cie_convergence() -------------------------------------------
% Example of the common information estimate (CIE) convergence. 
%
% The example is a decentralized data fusion scenario where two agents, 
% each equipped with a local sensor and a Kalman filter, in 2D measure and 
% track a common stationary target. Linear measurement models are used and
% the process noise is zero (can be changed). The agents communicate 
% dimension-reduced estimate derived using GEVO-KF. A decorrelation 
% procedure is applied to avoid double counting of information. KF is used 
% for fusion.
%
% The purpose with the example is to examine the convergence of the common
% information estimate, that is, to investigate if and how Gamma (G)
% converges to R1. The simulation is visualized using a plot of R1 vs Gamma
% for different k, and I vs inv(L)*Gamma*inv(L'), where R1=LL', for
% different k. With a process noise covariance Q=0 Gamma converges fast to 
% R1. The convergence depends on Q.
%
% The communication is according to the basic communication schedule.
% 
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- PARAMETERS ----------------------------------------------------------
nk = 8;
XT = zeros(2,nk);

pm.F = eye(2);
pm.Q = 0*eye(2);
pm.G0 = 100*eye(2);


% --- AGENTS --------------------------------------------------------------
agent = cell(2,1);

agent{1}.xhat = zeros(2,1); agent{1}.P = zeros(2);
agent{1}.ghat = zeros(2,1); agent{1}.G = zeros(2);
agent{1}.R = [16 4 ; 4 25];

agent{2}.xhat = zeros(2,1); agent{2}.P = zeros(2);
agent{2}.ghat = zeros(2,1); agent{2}.G = zeros(2);
agent{2}.R = [16 -4 ; -4 25];


% --- DATA RECORD ---------------------------------------------------------
data.R1 = cell(1,nk);
data.P = cell(1,nk);
data.G = cell(1,nk);
data.LGL = cell(1,nk);


% --- PLOT SETTINGS -------------------------------------------------------
close all
set_latex_interpreter;
f = figure(1);clf
f.Position = [100 100 1200 800];


% --- RUN SIMULATIONS -----------------------------------------------------
y = cell(2,1);
for k = 1:nk


    % --- 0. GENERATE MEAUREMENTS -----------------------------------------
    for i = 1:2; y{i} = generate_measurement(XT(:,k),agent{i}.R); end    


    % --- 1. INITIALIZE USING FIRST MEASUREMENT ---------------------------
    if k == 1
        for i = 1:2
            agent{i}.xhat = y{i};
            agent{i}.P = agent{i}.R;
            agent{i}.ghat = y{i};
            agent{i}.G = pm.G0;
        end


    % --- 2. PREDICT AND UDPATE LOCAL ESTIMATES ---------------------------    
    else
        for i = 1:2
            agent{i} = predict_estimates(agent{i},pm);
            agent{i} = update_estimates(agent{i},y{i});
        end
    end


    % --- 3. RUN GEVO AND SIMULATE COMMUNICATION --------------------------
    [itx,irx] = get_communication_topology(k);
    [agent{itx},tx_est] = run_gevo(agent{itx});


    % --- 4. FUSE RECEIVED ESTIMATE ---------------------------------------
    agent{irx} = fuse_received_estimate(agent{irx},tx_est);   


    % --- PLOT ------------------------------------------------------------
    LGL = plot_compare_ellipses(f,agent,k);


    % --- SAVE DATA -------------------------------------------------------
    data.R1{k} = agent{1}.P;
    data.P{k} = agent{2}.P;
    data.G{k} = agent{2}.G;
    data.LGL{k} = LGL;

end

end


% --- SENSOR AND FILTER FUNCTIONALITY -------------------------------------
function y = generate_measurement(xt,R)
    L = chol(R,'lower');
    y = xt + L*randn(2,1);
end

function agent = predict_estimates(agent,pm)
    F = pm.F; Q = pm.Q;
    agent.xhat = F*agent.xhat;
    agent.P = F*agent.P*F' + Q;
    agent.ghat = F*agent.ghat;
    agent.G = F*agent.G*F' + Q;   
end

function agent = update_estimates(agent,y)
    H = eye(2); R = agent.R; I = eye(2);
    ybar = y - H*agent.xhat;
    S = H*agent.P*H' + R;
    K = agent.P*H'/S;
    agent.xhat = agent.xhat + K*ybar;
    agent.P = (I - K*H)*agent.P;
end

function agent = fuse_received_estimate(agent,rx_est)
    [agent.xhat,agent.P] = fuse_kf(agent.xhat,agent.P,rx_est);
    [agent.ghat,agent.G] = fuse_kf(agent.ghat,agent.G,rx_est);
end

function [xhat,P] = fuse_kf(y1,R1,est)
    Psi = est.Psi; yPsi = est.yPsi; RPsi = est.RPsi;
    I1 = inv(R1); IPsi = inv(RPsi);
    P = inv(I1 + Psi'*IPsi*Psi);
    xhat = P*(I1*y1 + Psi'*IPsi*yPsi);
end


% --- GEVO FUNCTIONALITY --------------------------------------------------
function [agent,tx_est] = run_gevo(agent)
    y1 = agent.ghat; R1 = agent.G;
    [y2,R2] = subtract_estimate(agent.xhat,agent.P,y1,R1);
    Psi = gevo_kf(R1,R2,1);
    tx_est.Psi = Psi;
    tx_est.yPsi = Psi*y2;
    tx_est.RPsi = Psi*R2*Psi';
    [agent.ghat,agent.G] = fuse_kf(y1,R1,tx_est);
end

function Psi = gevo_kf(R1,R2,m)
    Q = R1*R1;
    S = R1+R2; 
    [X,D] = eig(Q,S,'qz');
    idx_vec = get_max_idx_vec(D,m);
    Psi = gram_schmidt_process(X(:,idx_vec))';
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

function [itx,irx] = get_communication_topology(k)
    if mod(k,2) == 0; itx = 2; irx = 1;
    else; itx = 1; irx = 2;
    end
end

function [c,C] = subtract_estimate(a,A,b,B)
    Ai = inv(A); Bi = inv(B);
    C = inv(Ai - Bi);
    c = C*(Ai*a - Bi*b);
end


% --- PLOT FUNCTIONALITY --------------------------------------------------
function c = get_colors()
    clr = get_thesis_colors;
    c.R1 = clr.darkyellow;
    c.G = clr.blue;
end

function varargout = plot_ellipse(x0,S,varargin)
    L = chol(S,'lower');
    N = 1000;
    a = linspace(0,2*pi,N);
    x = L*[cos(a) ; sin(a)];
    h = plot(x0(1)+x(1,:),x0(2)+x(2,:),varargin{:});  
    if nargout > 0; varargout{1} = h; end

    set_fontsize_all(14)
end

function LGL = plot_compare_ellipses(f,agent,k)
    x0 = [k;0];
    clr = get_colors; lw = 1.5;
    R1 = agent{1}.P; G = agent{2}.G;
    Li = inv(chol(R1,'lower')); LGL = Li*G*Li';

    figure(f)

    subplot(2,1,1);hold on
    h1 = plot_ellipse(10*x0,R1); h1.Color = clr.R1; h1.LineWidth = lw;
    hg = plot_ellipse(10*x0,G); hg.Color = clr.G; hg.LineWidth = lw;
    remove_ticks_and_ticklabels;     
    legend(gca,[h1 hg],'$R_1=LL^T$','$\Gamma$')
    axis equal; axis off

    subplot(2,1,2);hold on
    h1 = plot_ellipse(4.5*x0,eye(2)); h1.Color = clr.R1; h1.LineWidth = lw;
    hg = plot_ellipse(4.5*x0,LGL); hg.Color = clr.G; hg.LineWidth = lw;
    remove_ticks_and_ticklabels
    legend(gca,[h1 hg],'$I$','$L^{-1}\Gamma L^{-T}$')  
    axis equal; axis off

    set_fontsize_all(14)
end