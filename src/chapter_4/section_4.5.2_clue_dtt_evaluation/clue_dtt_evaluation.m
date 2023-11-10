function [agent,process] = clue_dtt_evaluation(varargin)
% --- clue_dtt_evaluation() -----------------------------------------------
% Section 4.5.2 Decentralized Target Tracking Using CLUE: Figure 4.10-4.11
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- PARAMETERS ----------------------------------------------------------
par.nagents = 3;
par.nk = 15;
par.M = 1000;
par.network_connectivity = 2;                                               % 1 = reduced, 2 = full
         
if nargin >= 2
    agent = varargin{1}; process = varargin{2}; 
    par.skip_precomp = 1;
else
    par.skip_precomp = 0;
end


% --- AGENTS --------------------------------------------------------------
if ~par.skip_precomp
    pref = 0.1;
    vref = 10;
    q = 2;
    process = get_process_model(pref,vref,q);

    agent = cell(par.nagents,1);
    for i = 1:par.nagents; agent{i} = get_agent_struct(par); end
    
    C0 = [100 0 ; 0 25]; T60 = [cosd(60) -sind(60) ; sind(60) cosd(60)];
    agent{1}.sensor.C = C0; 
    agent{2}.sensor.C = T60*C0*T60'; 
    agent{3}.sensor.C = T60'*C0*T60;

    for i = 1:par.nagents; agent{i}.sensor.L = chol(agent{i}.sensor.C,'lower'); end
end


% --- PRECOMPUTATIONS -----------------------------------------------------
fprintf('\n--- DTT USING CLUE ---\n\n')
fprintf('nagents = %d\n',par.nagents)
fprintf('nk = %d\n',par.nk)
fprintf('M = %d\n\n',par.M)

switch par.network_connectivity
    case 1; fprintf('Network connnectivity: reduced\n')
    case 2; fprintf('Network connnectivity: full\n')
end

[idx_tx,idx_rx] = get_basic_communication_schedule(par);

if ~par.skip_precomp
    fprintf('\n--- PRECOMPUTATIONS ---\n')

    for k = 1:par.nk
    
        fprintf('k = %d\n',k)
        
        if k == 1
            % INITIALIZE:
            for i = 1:par.nagents
                P0 = blkdiag(agent{i}.sensor.C,process.P0_vel);
                agent{i} = initiate_all_covariances(agent{i},P0);
            end
        else
            % LOCAL FILTERING:
            for i = 1:par.nagents
                agent{i} = time_update_all_covariances(agent{i},process,k);
                agent{i} = measurement_update_all_covariances(agent{i},k);
            end  
        end
    
        % COMMUNICATION:
        data_tx = agent{idx_tx{k}};
    
        % FUSION:
        for i = idx_rx{k}
            agent{i} = fusion_update_all_covariances(agent{i},data_tx,k);
        end
    
         % SAVE STUFF:
        for i = 1:par.nagents
            agent{i} = save_all_covariances(agent{i},k);
        end
    end
end


% --- ESTIMATION LOOP -----------------------------------------------------
rec = get_datarecord(par,process);

k0 = 1;
k25 = round(0.25*par.M);
k50 = round(0.5*par.M);
k75 = round(0.75*par.M);
k100 = par.M;

fprintf('\n--- ESTIMATION ---\n')
for m = 1:par.M

    switch m
        case k0; fprintf('0%%\n')
        case k25; fprintf('25%%\n')
        case k50; fprintf('50%%\n')
        case k75; fprintf('75%%\n')
    end

    % SIMULATE TARGET TRAJECTORY:
    xt = simulate_target_trajectory(par,process);

    % RUN ESTIMATION:
    for k = 1:par.nk
        
        if k == 1
            % INITIALIZE:
            for i = 1:par.nagents
                z = generate_target_measurement(agent{i},xt(:,k));
                agent{i} = initialize_local_estimates(agent{i},process,z);
            end
        else
            % LOCAL FILTERING:
            for i = 1:par.nagents
                z = generate_target_measurement(agent{i},xt(:,k));
                agent{i} = time_update_local_estimates(agent{i},process,k);
                agent{i} = measurement_update_local_estimates(agent{i},z,k);
            end  
        end

        % COMMUNICATION:
        data_tx = agent{idx_tx{k}};

        % FUSION:
        for i = idx_rx{k}
            agent{i} = fusion_update_local_estimates(agent{i},data_tx,k);
        end

        % SAVE STUFF:
        rec = save_local_estimates(rec,agent,xt,m,k);
    end

    if m == k100; fprintf('100%%\n\n'); end
end

% --- COMPUTE STATISTICS --------------------------------------------------
stats = compute_stats(rec,agent);


% --- PLOT ----------------------------------------------------------------
for i = 1:par.nagents; plot_stats(stats,i,par); end


end




% --- EVALUATION FUNCTIONS ------------------------------------------------
function stats = compute_stats(rec,agent)
    xt = rec.xt; 
    nagents = length(agent);
    M = length(xt);
    [nx,nk] = size(xt{1});
    a = zeros(1,nk);    
    stats = cell(nagents,1);

    f = get_method_fieldnames();
    for iagent = 1:nagents
        for i = 1:length(f)
            ag = agent{iagent}.(f{i}); 
            rec_ag = rec.agent{iagent}.(f{i});
            coin.(f{i}) = a; anees.(f{i}) = a; 
            rmt.(f{i}) = a; rmse.(f{i}) = a;
            S.(f{i}) = cell(1,nk);
            for k = 1:nk
                E = zeros(nx);    
                for m = 1:M
                    e = rec_ag.xhat{m}(:,k) - xt{m}(:,k);
                    E = E + e*e';
                end
                E = E/M;
                P = ag.P{k};
                S.(f{i}){k} = E;
                coin.(f{i})(k) = compute_coin(P,E);
                anees.(f{i})(k) = trace(E/P)/nx;
                rmt.(f{i})(k) = sqrt(trace(P));
                rmse.(f{i})(k) = sqrt(trace(E));
            end
        end
        stats{iagent}.coin = coin;
        stats{iagent}.anees = anees;
        stats{iagent}.rmt = rmt;
        stats{iagent}.rmse = rmse;
        stats{iagent}.S = S;
    end
end

function coin = compute_coin(P,E)
    L = chol(P,'lower'); 
    Li = inv(L);
    coin = max(eig(make_symmetric(Li*E*Li')));
end

function plot_stats(stats,iagent,par)
    coin = stats{iagent}.coin;    
    anees = stats{iagent}.anees;
    rmt = stats{iagent}.rmt;
    rmse = stats{iagent}.rmse;

    kvec = 1:par.nk;
    k_fu = get_time_of_fusion(par);
    idx = k_fu{iagent};

    clr = get_thesis_colors;
    lw = 2;

    figure(iagent);clf

    subplot(2,2,1);hold on
    hkf = plot(kvec(idx),coin.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),coin.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hici = plot(kvec(idx),coin.ici(idx),'-'); hici.Color = clr.ici; hici.LineWidth = lw;
    hle = plot(kvec(idx),coin.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    hro = plot(kvec(idx),coin.ro(idx),'-'); hro.Color = clr.ro; hro.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('COIN'); box on

    subplot(2,2,2);hold on
    hkf = plot(kvec(idx),anees.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),anees.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hici = plot(kvec(idx),anees.ici(idx),'-'); hici.Color = clr.ici; hici.LineWidth = lw;
    hle = plot(kvec(idx),anees.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    hro = plot(kvec(idx),anees.ro(idx),'-'); hro.Color = clr.ro; hro.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('ANEES'); box on

    subplot(2,2,3);hold on
    hkf = plot(kvec(idx),rmt.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),rmt.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hici = plot(kvec(idx),rmt.ici(idx),'-'); hici.Color = clr.ici; hici.LineWidth = lw;
    hle = plot(kvec(idx),rmt.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    hro = plot(kvec(idx),rmt.ro(idx),'-'); hro.Color = clr.ro; hro.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('RMT'); box on    

    subplot(2,2,4);hold on
    hkf = plot(kvec(idx),rmse.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),rmse.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hici = plot(kvec(idx),rmse.ici(idx),'-'); hici.Color = clr.ici; hici.LineWidth = lw;
    hle = plot(kvec(idx),rmse.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    hro = plot(kvec(idx),rmse.ro(idx),'-'); hro.Color = clr.ro; hro.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('RMSE'); box on

    legend(gca,[hkf hci hici hle hro],'KF','CI','ICI','LE','RO')

    set_fontsize_all(14)    

end





% --- TRACKING FUNCTIONS --------------------------------------------------
function z = generate_target_measurement(agent,xt)
    H = agent.sensor.H; L = agent.sensor.L; ny = size(H,1);
    z = H*xt + L*randn(ny,1);
end

function agent = initialize_local_estimates(agent,process,z)
    nvel = size(process.P0_vel,1);
    xhat0 = [z ; zeros(nvel,1)];
    f = get_method_fieldnames();
    for i = 1:length(f)
        agent.(f{i}).est.xhat = xhat0;
        agent.(f{i}).est.P = agent.(f{i}).P{1};
    end
end

function agent = time_update_local_estimates(agent,process,k)
    f = get_method_fieldnames();
    for i = 1:length(f)
        agent.(f{i}).est.xhat = process.F*agent.(f{i}).est.xhat;
        agent.(f{i}).est.P = agent.(f{i}).Ptu{k};
    end
end

function agent = measurement_update_local_estimates(agent,z,k)
    C = agent.sensor.C; H = agent.sensor.H;
    f = get_method_fieldnames();
    for i = 1:length(f)
        xhat = agent.(f{i}).est.xhat; P = agent.(f{i}).est.P; 
        K = P*H'/(H*P*H'+C);
        agent.(f{i}).est.xhat = xhat + K*(z-H*xhat);
        agent.(f{i}).est.P = agent.(f{i}).Pmu{k};
    end
end

function agent = fusion_update_local_estimates(agent,data_tx,k)
    f = get_method_fieldnames();
    for i = 1:length(f)
        y = [agent.(f{i}).est.xhat ; data_tx.(f{i}).est.xhat];
        K = agent.(f{i}).K{k}; P = agent.(f{i}).P{k};
        agent.(f{i}).est.xhat = K*y;
        agent.(f{i}).est.P = P;
    end
end

function rec = save_local_estimates(rec,agent,xt,m,k)
    f = get_method_fieldnames();
    for i = 1:length(agent)
        for j = 1:length(f)
            rec.agent{i}.(f{j}).xhat{m}(:,k) = agent{i}.(f{j}).est.xhat;
        end
    end
    rec.xt{m} = xt;
end




% --- PRECOMPUTATION FUNCTIONS --------------------------------------------
function agent = initiate_all_covariances(agent,P0)
    f = get_method_fieldnames();
    for i = 1:length(f)
        agent.(f{i}).est.P = P0;
        agent.(f{i}).Ptu{1} = P0;
        agent.(f{i}).Pmu{1} = P0;
    end
end

function agent = time_update_all_covariances(agent,process,k)
    F = process.F; Q = process.Q;
    f = get_method_fieldnames();
    for i = 1:length(f)
        P = agent.(f{i}).est.P;
        agent.(f{i}).est.P = F*P*F' + Q;
        agent.(f{i}).Ptu{k} = agent.(f{i}).est.P;
    end
end

function agent = measurement_update_all_covariances(agent,k)
    C = agent.sensor.C; H = agent.sensor.H;
    f = get_method_fieldnames();
    for i = 1:length(f)
        IP = inv(agent.(f{i}).est.P);
        agent.(f{i}).est.P = inv(IP + H'/C*H);
        agent.(f{i}).Pmu{k} = agent.(f{i}).est.P;
    end
end

function agent = fusion_update_all_covariances(agent,data_tx,k)
    f = get_method_fieldnames();
    for i = 1:length(f)
        R1 = agent.(f{i}).est.P; R2 = data_tx.(f{i}).est.P;
        switch lower(f{i})
            case 'kf'; [K,P] = kf_gain_and_covariance(R1,R2);                
            case 'ci'; [K,P] = ci_gain_and_covariance(R1,R2); 
            case 'ici'; [K,P] = ici_gain_and_covariance(R1,R2); 
            case 'le'; [K,P] = le_gain_and_covariance(R1,R2); 
            case 'ro'; [K,P] = ro_gain_and_covariance(R1,R2); 
            otherwise; error('unknown method...')
        end
        agent.(f{i}).K{k} = K;
        agent.(f{i}).P{k} = P;
        agent.(f{i}).est.P = P;
    end
end

function agent = save_all_covariances(agent,k)
    f = get_method_fieldnames();
    for i = 1:length(f)
        agent.(f{i}).P{k} = agent.(f{i}).est.P;
    end
end




% --- FUSION METHODS ------------------------------------------------------
function [K,P] = kf_gain_and_covariance(R1,R2)
    I1 = inv(R1); I2 = inv(R2);
    P = inv(I1 + I2);
    K = [P*I1 P*I2];
end

function [K,P] = ci_gain_and_covariance(R1,R2)
    I1 = inv(R1); I2 = inv(R2);
    J = @(w) trace(inv(w*I1 + (1-w)*I2));
    w = fminbnd(J,0,1);
    P = inv(w*I1 + (1-w)*I2);
    K = [w*P*I1 (1-w)*P*I2];
end

function [K,P] = ici_gain_and_covariance(R1,R2)
    I1 = inv(R1); I2 = inv(R2);
    J = @(w) trace(inv(I1 + I2 - inv(w*R1+(1-w)*R2)));
    w = fminbnd(J,0,1);
    IG = inv(w*R1+(1-w)*R2);
    P = inv(I1 + I2 - IG);
    K = [P*(I1-w*IG) P*(I2-(1-w)*IG)];
end

function [K,P] = le_gain_and_covariance(R1,R2)
    nx = size(R1,1);
    [U1,D1] = eig(R1); T1 = inv(sqrt(D1))*U1';
    [U2,D2] = eig(T1*R2*T1'); T2 = U2';
    T = T2*T1; Ti = inv(T);
    D = zeros(nx); K1 = zeros(nx);
    for i = 1:nx
        if D2(i,i) >= 1; K1(i,i) = 1; D(i,i) = 1;
        else; D(i,i) = D2(i,i);
        end
    end
    K2 = eye(nx) - K1;
    P = make_symmetric(Ti*D*Ti');
    K = [Ti*K1*T Ti*K2*T];
end

function [K,P] = ro_gain_and_covariance(R1,R2)
    %R1 = make_symmetric(R1); R2 = make_symmetric(R2);
    rhomax = 0.5;
    nx = size(R1,1);
    H = [eye(nx) ; eye(nx)];  
    P = sdpvar(nx);                                                                                                                      
    K = sdpvar(nx,2*nx);                                                                                                               
    R12 = sdpvar(nx,nx,'full');                                                   
    R = [R1 R12 ; R12' R2];   
    F = [K*H == eye(nx), uncertain(R12), [rhomax^2*R2 R12' ; R12 R1] >= 0];             
    F = [F, R >= 0, [P K*R ; R*K' R] >= 0];                                  
    J = trace(P);                                                           
    options = sdpsettings('solver','mosek','verbose',0,'debug',0); % ,'verbose',0,'debug',0
    optimize(F,J,options);  
    K = value(K);
    P = make_symmetric(value(P));
end




% --- BASIC FUNCTIONS -----------------------------------------------------
function xt = simulate_target_trajectory(par,process)
    nx = size(process.F,1);
    xt = zeros(nx,par.nk);
%     xt(:,1) = process.LP0*randn(nx,1);
    xt(:,1) = mvnrnd(zeros(1,nx),process.P0)';
    for k = 2:par.nk
        w = process.LQ*randn(nx,1);
        xt(:,k) = process.F*xt(:,k-1) + w;
    end
end

function A = get_agent_struct(par)
    c = cell(1,par.nk);
    s.K = c; s.P = c; s.Ptu = c; s.Pmu = c; s.est.xhat = []; s.est.P = [];
    f = get_method_fieldnames();
    for i = 1:length(f); A.(f{i}) = s; end
    A.sensor.H = [eye(2) zeros(2)]; A.sensor.R = [];
end

function process = get_process_model(pref,vref,q)
    Ts = 1; I2 = eye(2); O2 = zeros(2);
    F = [I2 Ts*I2 ; O2 I2]; 
    Q = q^2*[Ts^3*I2/3 Ts^2*I2/2 ; Ts^2*I2/2 Ts*I2];
    [U,D] = eig(Q);

    process.P0_pos = pref^2*eye(2);
    process.P0_vel = vref^2*eye(2);
    process.P0 = blkdiag(process.P0_pos,process.P0_vel); 
    process.LP0 = chol(process.P0,'lower');

    process.Ts = Ts;
    process.q = q;
    process.F = F;
    process.Q = Q;
    process.LQ = U*sqrt(D);
end

function [idx_tx,idx_rx] = get_basic_communication_schedule(par)
    T = get_commmunication_topology(par);
    idx_tx = cell(1,par.nk); idx_rx = idx_tx;
    for k = 1:par.nk
        idx_tx{k} = mod(k-1,par.nagents)+1;
        idx_rx{k} = find(T(idx_tx{k},:) == 1);
    end
end

function k_fu = get_time_of_fusion(par)
    [~,idx_rx] = get_basic_communication_schedule(par);
    k_fu = cell(par.nagents,1);
    for i = 1:par.nagents
        k_fu{i} = false(1,par.nk);
        for k = 1:par.nk
            if is_element_in_vector(i,idx_rx{k}); k_fu{i}(k) = true; end
        end
    end
end

function T = get_commmunication_topology(par)
    switch par.network_connectivity
        case 1
            T = zeros(par.nagents);
            T(par.nagents,1) = 1; 
            for i = 1:par.nagents-1; T(i,i+1) = 1; end
        case 2
            T = ones(par.nagents);
            for i = 1:par.nagents; T(i,i) = 0; end
    end
end

function fields = get_method_fieldnames()
    fields = {'kf' ; 'ci' ; 'ici' ; 'le' ; 'ro'};
end

function rec = get_datarecord(par,process)
    nx = size(process.F,1);
    c = cell(par.M,1); a = zeros(nx,par.nk);
    for m = 1:par.M; c{m} = a; end
    rec.agent = cell(par.nagents,1);
    rec.xt = c;
    f = get_method_fieldnames();
    for i = 1:length(f); s.(f{i}).xhat = c; end
    for i = 1:par.nagents; rec.agent{i} = s; end
end

