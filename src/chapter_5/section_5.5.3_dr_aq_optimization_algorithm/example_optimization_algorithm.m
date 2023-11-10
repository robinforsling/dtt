function example_optimization_algorithm
% --- example_optimization_algorithm() ------------------------------------
% Sedction 5.5.3 Preserving Correction Assignment - Example Optimization
% Algorithm: Figure 5.19
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- PARAMETERS ---
M = 1000;
n = 4;
p = 2;

scenID = 4;


% --- DEFINE PROBLEM ------------------------------------------------------
[Y,S] = load_scenario(scenID,n,p);
p = length(Y); n = length(Y{1});


% --- INITIATE ------------------------------------------------------------
U = cell(p,1); Lambda = U;
for i = 1:p
    [X,E] = eig(Y{i},S{i});
    idx = get_max_idx_vec(E,1);
    Lambda{i} = E(idx,idx);
    u = X(:,idx); U{i} = u/norm(u);
end

z0 = zeros(n,1);
for i = 1:p; z0 = z0 + U{i}/sqrt(Lambda{i}); end % 
%z0 = randn(n,1);
z0 = z0/norm(z0);


% --- OPTIMIZATION SETTINGS -----------------------------------------------
opt.K = 25;
opt.z0 = z0;
opt.alpha = 0.001;
opt.alpha_low = 0.01;
opt.alpha_high = 1;
opt.alpha_rnd = 0.0;


% --- SIMULATE ------------------------------------------------------------

% LINESEARCH:
[Psi_ls,J_ls,Jmin_ls] = optimize_psi_linesearch(opt,U,Y,S);

% FIXED SMALL STEP SIZE:
opt.alpha = 0.01;
[Psi_small,J_small,Jmin_small] = optimize_psi_fixed_step_size(opt,U,Y,S);

% FIXED LARGE STEP SIZE:
opt.alpha = 1;
[Psi_large,J_large,Jmin_large] = optimize_psi_fixed_step_size(opt,U,Y,S);

Z = randn(n,M);
J_samp = zeros(p,M);
for k = 1:M
    z = Z(:,k); z = z/norm(z);
    for i = 1:p
        J_samp(i,k) = evaluate_q(z,Y{i},S{i});
    end
end


% --- PLOT ----------------------------------------------------------------
clr = get_thesis_colors;
clrsamp = clr.yellow; 
ms = 2;

figure(2);clf;

subplot(1,2,1);hold on
hs = plot(J_samp(1,:),J_samp(2,:),'o'); hs.Color = clrsamp; hs.MarkerSize = ms;
hsm = plot(J_small(1,:),J_small(2,:),'x'); hsm.Color = clr.sm;
hsmf = plot(J_small(1,end),J_small(2,end),'o'); hsmf.Color = clr.sm;
hlg = plot(J_large(1,:),J_large(2,:),'+'); hlg.Color = clr.lg;
hlgf = plot(J_large(1,end),J_large(2,end),'o'); hlgf.Color = clr.lg;
hls = plot(J_ls(1,:),J_ls(2,:),'*-'); hls.Color = clr.ls;
hlsf = plot(J_ls(1,end),J_ls(2,end),'o'); hlsf.Color = clr.ls;
legend(gca,[hsm hlg hls],'Small step size','Large step size','Adaptive step size')
axis equal

subplot(1,2,2);hold on
hsm = plot(Jmin_small,'-'); hsm.Color = clr.sm;
hlg = plot(Jmin_large,'-'); hlg.Color = clr.lg;
hls = plot(Jmin_ls,'-'); hls.Color = clr.ls;
legend(gca,[hsm hlg hls],'Small step size','Large step size','Adaptive step size')

set_fontsize_all(14)

end



% --- MISC FUNCTIONS ------------------------------------------------------
function [Psi,J,Jmin] = optimize_psi_linesearch(opt,U,Y,S)

    p = length(U);
    n = length(Y{1});

    V = eye(n);
    
    % Optimizer parameters:
    c = 1;
    K = opt.K;
    z = opt.z0;
    alpha_rnd = opt.alpha_rnd;
    alpha_low = opt.alpha_low;
    alpha_high = opt.alpha_high;

    % Optimize:
    Jmin = zeros(1,K); idxmin = Jmin; 
    J = zeros(p,K);
    for k = 1:K
    
        % Compute q:
        q = zeros(1,p); Duq = q;
        for i = 1:p; q(i) = evaluate_q(z,Y{i},S{i}); end
        J(:,k) = q';
    
        idx = get_min_idx_vec(q,1);
        Jmin(k) = q(idx); idxmin(k) = idx;
    
        % Directional derivatives:
        u = U{idx};
        %u = V(:,ceil(n*rand));
        for i = 1:p; Duq(i) = directional_derivative(z,u,q(i),Y{i},S{i}); end

        % Select step size:
        a = zeros(1,p); a(idx) = Inf;
        for i = get_index_vector(p,idx)
            a(i) = (q(i)-q(idx))/(Duq(idx)-Duq(i));
        end
        a(idx) = [];
    
        if Duq(idx) >= 0; alpha = min([a(a>0)]);
        else; alpha = max([a(a<0)]);
        end
    
        if isempty(alpha); alpha = alpha_low*sign(Duq(idx)); end
        if abs(alpha) < alpha_low; alpha = alpha_low*sign(Duq(idx)); end
        if abs(alpha) > alpha_high; alpha = alpha_high*sign(Duq(idx)); end
    
        u_rnd = U{ceil(p*rand)};
    
        % Update z:
        z = z + c*alpha*u + alpha_rnd*u_rnd;
        z = z/norm(z);
    end

    Psi = z';
end

function [Psi,J,Jmin] = optimize_psi_fixed_step_size(opt,U,Y,S)
    p = length(U);
    
    % Optimizer parameters:
    c = 1;
    K = opt.K;
    z = opt.z0;
    alpha_rnd = opt.alpha_rnd;

    % Optimize:
    Jmin = zeros(1,K); idxmin = Jmin; 
    J = zeros(p,K);
    for k = 1:K
    
        % Compute q:
        q = zeros(1,p); 
        for i = 1:p; q(i) = evaluate_q(z,Y{i},S{i}); end
        J(:,k) = q';
    
        idx = get_min_idx_vec(q,1);
        Jmin(k) = q(idx); idxmin(k) = idx;

        u = U{idx};
        u_rnd = U{ceil(p*rand)};
        
        Duq = directional_derivative(z,u,q(idx),Y{idx},S{idx});
        alpha = opt.alpha*sign(Duq); 
    
        % Update z:
        z = z + c*alpha*u + alpha_rnd*u_rnd;
        z = z/norm(z);
    end

    Psi = z';
end

function q = evaluate_q(z,Y,S)
    q = (z'*Y*z)/(z'*S*z);
end

function Du = directional_derivative(z,u,q,Y,S)
    dqdx = q_gradient(z,q,Y,S);
    Du = u'*dqdx;
end

function dqdx = q_gradient(z,q,Y,S)
    dqdx = 2*((Y-q*S)/(z'*S*z))*z;
end

function idx = get_index_vector(n,iremove)
    idx = 1:n;
    idx(iremove) = [];
end

function [Y,S] = load_scenario(scenID,n,p)

    switch scenID
    
        % CASE 1
        case 1
            if p > n
                a = 100*eye(n,1);
                for i = 1:p
                    R2{i} = get_random_covariance(n);
                    y = get_rnd_orthogonal_mat(n)*a;
                    Y{i} = y*y';
                    S{i} = 2*R2{i};
                end
            
            else
                a = eye(n,p);
                for i = 1:p
                    R2{i} = get_random_covariance(n);
                    y = 20*a(:,i);
                    Y{i} = y*y';
                    S{i} = 2*R2{i};
                end
            end
        
        case 3
            Y = cell(2,1); S = Y;
            n = 3;

            y1 = [50;40;-10];
            y2 = [60;120;10];

            Y{1} = y1*y1';
            Y{2} = y2*y2';

            W = diag([3 2 1]);
            S{1} = W*get_random_covariance(n)*W; 
            S{2} = W*get_random_covariance(n)*W; 

        case 4
            Y = cell(2,1); S = Y;
            n = 4;

            y1 = [50;40;-10;30];
            y2 = [60;120;-30;10];

            Y{1} = y1*y1';
            Y{2} = y2*y2';

            W = diag([5 5 1 1]);
            S{1} = W*get_random_covariance(n)*W; 
            S{2} = W*get_random_covariance(n)*W; 

        case 6
            Y = cell(2,1); S = Y;
            n = 6;

            y1 = [70;20;-10;60;0;0];
            y2 = [60;120;-30;10;0;0];

            Y{1} = y1*y1';
            Y{2} = y2*y2';

            W = diag([5 5 3 3 1 1]);
            S{1} = W*get_random_covariance(n)*W; 
            S{2} = W*get_random_covariance(n)*W; 
    end

end