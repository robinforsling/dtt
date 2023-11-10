function data = dr_aq_evaluation(varargin)
% --- dr_aq_evaluation() --------------------------------------------------
% Section 5.5.4 Numerical Evaluation - DR for Association Quality:
% Figure 5.20-5.21
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- PARAMETERS (DEFAULT) ------------------------------------------------
scenID = 2;
M = 10;
c_vec = 0.1:0.1:5;
nc = length(c_vec);

COVARIANCE_SAMPLING = 2;


% --- READ INPUTS ---------------------------------------------------------
if nargin > 0; scenID = varargin{1}; end
if nargin > 1; M = varargin{2}; end
if nargin > 2; c_vec = varargin{3}; nc = length(c_vec); end


% --- SCENARIO ------------------------------------------------------------
scen = load_scenario(scenID);
N = scen.N;
n = scen.n;


% --- AGENTS --------------------------------------------------------------
W = diag([3 3 2 2 1 1]);
A = cell(N,1); B = A;
switch COVARIANCE_SAMPLING
    case 0
        A0 = W*get_random_covariance(n)*W; 
        B0 = W*get_random_covariance(n)*W; 
        for i = 1:N
            A{i} = A0;
            B{i} = B0;
        end
    case 1
        for i = 1:N
            A{i} = W*get_random_covariance(n)*W;
            B{i} = W*get_random_covariance(n)*W;
        end
    case 2
        A0 = W*get_random_covariance(n)*W;
        B0 = W*get_random_covariance(n)*W;
        a0 = norm(A0); b0 = norm(B0); r0 = (a0+b0)/2;
        A0 = (r0/a0)*A0; B0 = (r0/b0)*B0;
        f = 0.75;
        for i = 1:N
            A{i} = A0 + f*W*get_random_covariance(n)*W;
            B{i} = B0 + f*W*get_random_covariance(n)*W;
        end
end


% --- MAIN LOOP -----------------------------------------------------------
CA_FULL = zeros(M,nc); CA_RED_F = CA_FULL; CA_RED_A = CA_FULL;
IA_FULL = zeros(M,nc); IA_RED_F = IA_FULL; IA_RED_A = IA_FULL;

for ic = 1:nc

    c = c_vec(ic);
    C = blkdiag(c*eye(2),eye(n-2));
    X = scen.X;   

    
    % --- COVARIANCES ---
    R1 = cell(N,1); R2 = cell(N,2);
    L1 = cell(N,1); L2 = cell(N,1);
    for i = 1:N
        switch COVARIANCE_SAMPLING
            case 0
                R1{i} = C*A*C; 
                R2{i} = C*B*C; 
            case {1,2}
                R1{i} = C*A{i}*C; 
                R2{i} = C*B{i}*C; 
        end
        L1{i} = chol(R1{i},'lower'); 
        L2{i} = chol(R2{i},'lower');
    end
    
    
    % --- FUSION OPTIMAL PSI ---
    PSI_F = cell(N,1); 
    for i = 1:N; PSI_F{i} = gevo_kf(R1{i},R2{i}); end
    
    
    % --- MC LOOP ---
    for k = 1:M
    
        Y1 = cell(N,1); Y2 = cell(N,1);    
        for i = 1:N
            Y1{i} = X{i} + L1{i}*randn(n,1);
            Y2{i} = X{i} + L2{i}*randn(n,1);
        end


        % --- ASSOCIATION BASED PSI ---
        PSI_A = cell(N,1); 
        for i = 1:N; PSI_A{i} = association_psi_linesearch(i,Y2,R2); end


        % --- SOLVE ASSIGNMENT PROBLEM ---
        AFULL = full_assignment_matrix(Y1,R1,Y2,R2);
        ARED_F = reduced_assignment_matrix(Y1,R1,Y2,R2,PSI_F);
        ARED_A = reduced_assignment_matrix(Y1,R1,Y2,R2,PSI_A);
        
        CFULL = solve_assignment_problem(AFULL);
        CRED_F = solve_assignment_problem(ARED_F);
        CRED_A = solve_assignment_problem(ARED_A);

        [CA_FULL(k,ic),IA_FULL(k,ic)] = evaluate_assignment_results(CFULL);
        [CA_RED_F(k,ic),IA_RED_F(k,ic)] = evaluate_assignment_results(CRED_F);
        [CA_RED_A(k,ic),IA_RED_A(k,ic)] = evaluate_assignment_results(CRED_A);
    end
end

R1 = cell(N,1); R2 = cell(N,2);
for i = 1:N
    switch COVARIANCE_SAMPLING
        case 0
            R1{i} = A; R2{i} = B; 
        case {1,2}
            R1{i} = A{i}; R2{i} = B{i}; 
    end
end

% --- PACK OUTPUT DATA ----------------------------------------------------
data.scen = scen;

data.sim.M = M;
data.sim.nc = nc;
data.sim.c_vec = c_vec;

data.est.R1 = R1;
data.est.R2 = R2;

data.full = compute_stats(CA_FULL,IA_FULL);
data.red_fus = compute_stats(CA_RED_F,IA_RED_F);
data.red_asso = compute_stats(CA_RED_A,IA_RED_A);


% --- VISUALIZE RESULTS --------------------------------------------------
data.prob = visualize_results(data);


end




% --- PSI SELECTION FUNCTIONS ---------------------------------------------
function Psi = association_psi_linesearch(itgt,Y2,R2)

    N = length(Y2); n = length(Y2{1}); p = N-1;
    y0 = Y2{itgt};
    idx = get_index_vector(N,itgt);
    Y2 = Y2(idx);
    R2 = R2(idx);

    % Initialize:
    S = cell(p,1);
    Y = cell(p,1);
    Lambda = zeros(1,p); U = cell(1,p);
    for i = 1:p
        y = Y2{i}-y0;
        Y{i} = y*y';
        S{i} = 2*R2{i};
        [X,E] = eig(Y{i},S{i});
        idx = get_max_idx_vec(E,1);
        Lambda(i) = E(idx,idx);
        U{i} = X(:,idx)/norm(X(:,idx));
    end
    
    z0 = zeros(n,1);
    for i = 1:p; z0 = z0 + U{i}/Lambda(i); end
    z = z0/norm(z0);

    K = 100;
    alpha_rnd = 0.0;
    alpha_low = 0.001;
    alpha_high = 1;
    qmin = zeros(1,K); idxmin = qmin; 
    for k = 1:K
    
        % Compute q:
        q = zeros(1,p); Duq = q;
        for i = 1:p; q(i) = evaluate_q(z,Y{i},S{i}); end
    
        idx = get_min_idx_vec(q,1);
        qmin(k) = q(idx); idxmin(k) = idx;
    
        % Directional derivatives:
        u = U{idx};
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
        z = z + alpha*u + alpha_rnd*u_rnd;
        z = z/norm(z);
    end

    Psi = z';
end

function Psi = gevo_kf(R1,R2)
    Q = R1*R1;
    S = R1+R2;
    [X,D] = eig(Q,S);
    idx = get_max_idx_vec(D,1);
    Psi = X(:,idx)';
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




% --- ASSIGNMENT FUNCIONS -------------------------------------------------
function C = solve_assignment_problem(A)
    max_cost = 1e10;
    C = matchpairs(A,max_cost);
    [~,idx] = sort(C(:,1));
    C = C(idx,:); 
end

function A = full_assignment_matrix(Y1,R1,Y2,R2)
    N = length(Y1);
    A = zeros(N);
    for i = 1:N
        for j = 1:N
            A(i,j) = md_full(Y1{i},R1{i},Y2{i},R2{i});
        end
    end
end

function A = reduced_assignment_matrix(Y1,R1,Y2,R2,PSI)
    N = length(Y1);
    A = zeros(N);
    for i = 1:N
        for j = 1:N
            A(i,j) = md_reduced(Y1{i},R1{i},Y2{j},R2{j},PSI{j});
        end
    end
end

function d2 = md_full(y1,R1,y2,R2)
    dy = y1-y2; 
    S = R1+R2;
    d2 = dy'/S*dy;
end

function r2 = md_reduced(y1,R1,y2,R2,Psi)
    dy = Psi*(y1-y2); 
    S = Psi*(R1+R2)*Psi';
    r2 = dy'/S*dy;
end




% --- EVALUATION FUNCTIONS ------------------------------------------------
function [nca,nia] = evaluate_assignment_results(C)
    N = size(C,1);
    nca = 0;
    for i = 1:N 
        if C(i,1) == C(i,2); nca = nca + 1; end
    end
    nia = N-nca;
end

function s = compute_stats(CA,IA)
    [M,nc] = size(CA);
    s.ca.mean = zeros(1,nc); s.ca.std = zeros(1,nc);
    s.ia.mean = zeros(1,nc); s.ia.std = zeros(1,nc);
    for ic = 1:nc
        s.ca.mean(ic) = mean(CA(:,ic));
        s.ca.std(ic) = std(CA(:,ic));
        s.ia.mean(ic) = mean(IA(:,ic));
        s.ia.std(ic) = std(IA(:,ic));
    end
end



% --- SCENARIO ------------------------------------------------------------
function scen = load_scenario(scen_ID)

    scen.ID = scen_ID;
    scen.N = [];
    scen.n = [];
    scen.X = [];
    
    switch scen_ID
    
        case 1
            N = 6;
            n = 6;
            X = cell(N,1); for i = 1:N; X{i} = zeros(n,1); end
            X{1}(1:4) = [0;0;40;-50];
            X{2}(1:4) = [10;100;50;-40];
            X{3}(1:4) = [50;50;-10;60];
            X{4}(1:4) = [60;120;-30;10];
            X{5}(1:4) = [90;20;-30;-10];
            X{6}(1:4) = [110;40;-20;50];
    
        case 2
            N = 10;
            n = 6;
            X = cell(N,1); for i = 1:N; X{i} = zeros(n,1); end
            X{1}(1:4) = [0;0;30;-50];
            X{2}(1:4) = [10;100;50;-30];
            X{3}(1:4) = [50;50;-10;60];
            X{4}(1:4) = [60;120;-30;10];
            X{5}(1:4) = [90;80;-30;-10];
            X{6}(1:4) = [110;40;-20;50];
            X{7}(1:4) = [30;-10;30;-30];
            X{8}(1:4) = [80;-10;20;40];
            X{9}(1:4) = [-20;80;30;-10];
            X{10}(1:4) = [-30;30;0;50];
     
        otherwise 
    end
    
    scen.N = N;
    scen.n = n;
    scen.X = X;
end




% --- MISC FUNCTIONS ------------------------------------------------------
function idx = get_index_vector(n,iremove)
    idx = 1:n;
    idx(iremove) = [];
end

