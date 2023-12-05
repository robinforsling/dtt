function [rec,agent,stats] = gevo_dtt_evaluation
% --- gevo_dtt_evaluation() -----------------------------------------------
% DTT scenario in the numerical evaluation of GEVO in the paper:
% "Decentralized State Estimation In A Dimension-Reduced Linear Regression"
%
% 2023-11-16 Robin Forsling

addpath(genpath(pwd))

set_latex_interpreter;


% --- PARAMETERS ----------------------------------------------------------
par.m = 2;
par.nagents = 3;
par.nk = 15;
par.M = 10;
                                                        

% --- AGENTS --------------------------------------------------------------
pref = 0.1;
vref = 10;
q = 2;
process = get_process_model(pref,vref,q);

agent = cell(par.nagents,1);
for i = 1:par.nagents; agent{i} = get_agent_struct(par); end

C0 = [100 0 ; 0 10]; T60 = [cosd(60) -sind(60) ; sind(60) cosd(60)];
agent{1}.sensor.C = C0; 
agent{2}.sensor.C = T60*C0*T60'; 
agent{3}.sensor.C = T60'*C0*T60;

for i = 1:par.nagents; agent{i}.sensor.L = chol(agent{i}.sensor.C,'lower'); end


% --- ESTIMATION LOOP -----------------------------------------------------
fprintf('\n--- GEVO EVALUATION ---\n\n')
fprintf('m = %d\n',par.m)
fprintf('nagents = %d\n',par.nagents)
fprintf('nk = %d\n',par.nk)
fprintf('M = %d\n\n',par.M)
fprintf('Network connnectivity: reduced\n')


[idx_tx,idx_rx] = get_basic_communication_schedule(par);

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
                agent{i} = time_update_local_estimates(agent{i},process);
                agent{i} = measurement_update_local_estimates(agent{i},z);
            end  
        end

        % GEVO & COMMUNICATION:
        agent_tx = agent{idx_tx{k}}; agent_rx = agent{idx_rx{k}};
        data_tx = manage_communication(par,agent_tx,agent_rx);

        % FUSION:
        for i = idx_rx{k}
            agent{i} = fusion_update_local_estimates(agent{i},data_tx);
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

    % DR:
    for iagent = 1:nagents
        for i = 1:length(f)
            rec_ag = rec.agent{iagent}.dr.(f{i});
            coin.(f{i}) = a; anees.(f{i}) = a; rmt.(f{i}) = a; 
            S.(f{i}) = cell(1,nk);
            for k = 1:nk
                E = zeros(nx);    
                for m = 1:M
                    e = rec_ag.xhat{m}(:,k) - xt{m}(:,k);
                    E = E + e*e';
                end
                E = E/M;
                P = rec_ag.P(:,:,k);
                S.(f{i}){k} = E;
                coin.(f{i})(k) = compute_coin(P,E);
                anees.(f{i})(k) = trace(E/P)/nx;
                rmt.(f{i})(k) = sqrt(trace(P));
            end
        end
        stats{iagent}.dr.coin = coin;
        stats{iagent}.dr.anees = anees;
        stats{iagent}.dr.rmt = rmt;
        stats{iagent}.dr.S = S;
    end

    % FULL:
    for iagent = 1:nagents
        for i = 1:length(f)
            rec_ag = rec.agent{iagent}.full.(f{i});
            coin.(f{i}) = a; anees.(f{i}) = a; rmt.(f{i}) = a; 
            S.(f{i}) = cell(1,nk);
            for k = 1:nk
                E = zeros(nx);    
                for m = 1:M
                    e = rec_ag.xhat{m}(:,k) - xt{m}(:,k);
                    E = E + e*e';
                end
                E = E/M;
                P = rec_ag.P(:,:,k);
                S.(f{i}){k} = E;
                coin.(f{i})(k) = compute_coin(P,E);
                anees.(f{i})(k) = trace(E/P)/nx;
                rmt.(f{i})(k) = sqrt(trace(P));
            end
        end
        stats{iagent}.full.coin = coin;
        stats{iagent}.full.anees = anees;
        stats{iagent}.full.rmt = rmt;
        stats{iagent}.full.S = S;
    end

    % REF:
    for iagent = 1:nagents
        rec_ag = rec.agent{iagent}.ref.dca;
        coin = a; anees = a; rmt = a; S = cell(1,nk);
        for k = 1:nk
            E = zeros(nx);    
            for m = 1:M
                e = rec_ag.xhat{m}(:,k) - xt{m}(:,k);
                E = E + e*e';
            end
            E = E/M;
            P = rec_ag.P(:,:,k);
            S{k} = E;
            coin(k) = compute_coin(P,E);
            anees(k) = trace(E/P)/nx;
            rmt(k) = sqrt(trace(P));
        end
        stats{iagent}.ref.dca.coin = coin;
        stats{iagent}.ref.dca.anees = anees;
        stats{iagent}.ref.dca.rmt = rmt;
        stats{iagent}.ref.dca.S = S;
    end
    for iagent = 1:nagents
        rec_ag = rec.agent{iagent}.ref.full;
        coin = a; anees = a; rmt = a; S = cell(1,nk);
        for k = 1:nk
            E = zeros(nx);    
            for m = 1:M
                e = rec_ag.xhat{m}(:,k) - xt{m}(:,k);
                E = E + e*e';
            end
            E = E/M;
            P = rec_ag.P(:,:,k);
            S{k} = E;
            coin(k) = compute_coin(P,E);
            anees(k) = trace(E/P)/nx;
            rmt(k) = sqrt(trace(P));
        end
        stats{iagent}.ref.full.coin = coin;
        stats{iagent}.ref.full.anees = anees;
        stats{iagent}.ref.full.rmt = rmt;
        stats{iagent}.ref.full.S = S;
    end

    % RATIO:
    for iagent = 1:nagents
        for i = 1:length(f)
            stats{iagent}.ratio.rmt.(f{i}) = stats{iagent}.dr.rmt.(f{i}) ./ stats{iagent}.full.rmt.(f{i});
        end
        stats{iagent}.ratio.rmt.ref = stats{iagent}.ref.dca.rmt ./ stats{iagent}.ref.full.rmt;
    end
end

function coin = compute_coin(P,E)
    L = chol(P,'lower'); 
    Li = inv(L);
    coin = max(eig(make_symmetric(Li*E*Li')));
end

function plot_stats(stats,iagent,par)
    coin = stats{iagent}.dr.coin; anees = stats{iagent}.dr.anees; rmtr = stats{iagent}.ratio.rmt;
    coin_ref = stats{iagent}.ref.dca.coin; anees_ref = stats{iagent}.ref.dca.anees; 

    kvec = 1:par.nk;
    k_fu = get_time_of_fusion(par);
    idx = k_fu{iagent};

    clr = get_paper_colors;
    lw = 2;

    figure(iagent);clf

    subplot(3,1,1);hold on
    hkf = plot(kvec(idx),coin.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),coin.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hle = plot(kvec(idx),coin.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    href = plot(kvec(idx),coin_ref(idx),'k--'); href.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('COIN'); box on

    subplot(3,1,2);hold on
    hkf = plot(kvec(idx),anees.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),anees.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hle = plot(kvec(idx),anees.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    href = plot(kvec(idx),anees_ref(idx),'k--'); href.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('ANEES'); box on

    subplot(3,1,3);hold on
    hkf = plot(kvec(idx),rmtr.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),rmtr.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hle = plot(kvec(idx),rmtr.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    href = plot(kvec(idx),rmtr.ref(idx),'k--'); href.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('RMTR'); box on    

    legend(gca,[hkf hci hle href],'KF','CI','LE','DCA')

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
    P0 = blkdiag(agent.sensor.C,process.P0_vel);
    f = get_method_fieldnames();
    for i = 1:length(f)
        agent.dr.(f{i}).est.xhat = xhat0;
        agent.dr.(f{i}).est.P = P0;
        agent.full.(f{i}).est.xhat = xhat0;
        agent.full.(f{i}).est.P = P0;
    end
    % REF:
    agent.ref.dca.est.xhat = xhat0;
    agent.ref.dca.est.P = P0;
    agent.ref.full.est.xhat = xhat0;
    agent.ref.full.est.P = P0;
end

function agent = time_update_local_estimates(agent,process)
    F = process.F; Q = process.Q;
    f = get_method_fieldnames();
    for i = 1:length(f)
        agent.dr.(f{i}).est.xhat = F*agent.dr.(f{i}).est.xhat;
        agent.dr.(f{i}).est.P = F*agent.dr.(f{i}).est.P*F' + Q;
        agent.full.(f{i}).est.xhat = F*agent.full.(f{i}).est.xhat;
        agent.full.(f{i}).est.P = F*agent.full.(f{i}).est.P*F' + Q;
    end
    % REF:
    agent.ref.dca.est.xhat = F*agent.ref.dca.est.xhat;
    agent.ref.dca.est.P = F*agent.ref.dca.est.P*F' + Q;

    agent.ref.full.est.xhat = F*agent.ref.full.est.xhat;
    agent.ref.full.est.P = F*agent.ref.full.est.P*F' + Q;
end

function agent = measurement_update_local_estimates(agent,z)
    C = agent.sensor.C; H = agent.sensor.H; 
    In = eye(size(H,2));
    f = get_method_fieldnames();
    for i = 1:length(f)
        % DR:
        xhat = agent.dr.(f{i}).est.xhat; P = agent.dr.(f{i}).est.P; 
        K = P*H'/(H*P*H'+C);
        agent.dr.(f{i}).est.xhat = xhat + K*(z-H*xhat);
        agent.dr.(f{i}).est.P = (In-K*H)*P;
        % FULL:
        xhat = agent.full.(f{i}).est.xhat; P = agent.full.(f{i}).est.P; 
        K = P*H'/(H*P*H'+C);
        agent.full.(f{i}).est.xhat = xhat + K*(z-H*xhat);
        agent.full.(f{i}).est.P = (In-K*H)*P;
    end
    % REF:
    xhat = agent.ref.dca.est.xhat; P = agent.ref.dca.est.P; 
    K = P*H'/(H*P*H'+C);
    agent.ref.dca.est.xhat = xhat + K*(z-H*xhat);
    agent.ref.dca.est.P = (In-K*H)*P;

    xhat = agent.ref.full.est.xhat; P = agent.ref.full.est.P; 
    K = P*H'/(H*P*H'+C);
    agent.ref.full.est.xhat = xhat + K*(z-H*xhat);
    agent.ref.full.est.P = (In-K*H)*P;
end

function agent = fusion_update_local_estimates(agent,data_tx)
    f = get_method_fieldnames();
    for i = 1:length(f)
        % DR:
        y1 = agent.dr.(f{i}).est.xhat; R1 = agent.dr.(f{i}).est.P;
        yPsi = data_tx.dr.(f{i}).est.xhat;
        RPsi = data_tx.dr.(f{i}).est.P;
        Psi = data_tx.dr.(f{i}).est.H;
        switch f{i}
            case 'kf'; [xhat,P] = kf_fusion(y1,R1,yPsi,RPsi,Psi);
            case 'ci'; [xhat,P] = ci_fusion(y1,R1,yPsi,RPsi,Psi);
            case 'le'; [xhat,P] = le_fusion(y1,R1,yPsi,RPsi,Psi);
        end
        agent.dr.(f{i}).est.xhat = xhat; agent.dr.(f{i}).est.P = P;
        % FULL:
        y1 = agent.full.(f{i}).est.xhat; R1 = agent.full.(f{i}).est.P;
        y2 = data_tx.full.(f{i}).est.xhat; R2 = data_tx.full.(f{i}).est.P;
        H2 = eye(length(y2));
        switch f{i}
            case 'kf'; [xhat,P] = kf_fusion(y1,R1,y2,R2,H2);
            case 'ci'; [xhat,P] = ci_fusion(y1,R1,y2,R2,H2);
            case 'le'; [xhat,P] = le_fusion(y1,R1,y2,R2,H2);
        end
        agent.full.(f{i}).est.xhat = xhat; agent.full.(f{i}).est.P = P;
    end
    % REF:
    y1 = agent.ref.dca.est.xhat; R1 = agent.ref.dca.est.P;
    y2 = data_tx.ref.dca.est.xhat; R2 = data_tx.ref.dca.est.P; H2 = eye(length(y2));
    [xhat,P] = ci_fusion(y1,R1,y2,R2,H2);
    agent.ref.dca.est.xhat = xhat; agent.ref.dca.est.P = P;

    y1 = agent.ref.full.est.xhat; R1 = agent.ref.full.est.P;
    y2 = data_tx.ref.full.est.xhat; R2 = data_tx.ref.full.est.P; H2 = eye(length(y2));
    [xhat,P] = ci_fusion(y1,R1,y2,R2,H2);
    agent.ref.full.est.xhat = xhat; agent.ref.full.est.P = P;
end

function data_tx = manage_communication(par,agent_tx,agent_rx)
    m = par.m;
    data_tx = agent_tx;

    % DR:
    f = get_method_fieldnames();
    for i = 1:length(f)
        R1 = agent_rx.dr.(f{i}).est.P;
        y2 = agent_tx.dr.(f{i}).est.xhat;
        R2 = agent_tx.dr.(f{i}).est.P;
        
        switch lower(f{i})
            case 'kf'; Psi = gevokf_dr(R1,R2,m);
            case 'ci'; Psi = gevoci_dr(R1,R2,m);
            case 'le'; Psi = gevole_dr(R1,R2,m); 
        end

        data_tx.dr.(f{i}).est.xhat = Psi*y2;
        data_tx.dr.(f{i}).est.P = Psi*R2*Psi';
        data_tx.dr.(f{i}).est.H = Psi;
    end

    % REF:
    data_tx.ref.dca.est.P = dca_eig(agent_tx.ref.dca.est.P);
end

function rec = save_local_estimates(rec,agent,xt,m,k)
    f = get_method_fieldnames();
    for i = 1:length(agent)
        for j = 1:length(f)
            rec.agent{i}.dr.(f{j}).xhat{m}(:,k) = agent{i}.dr.(f{j}).est.xhat;
            rec.agent{i}.dr.(f{j}).P(:,:,k) = agent{i}.dr.(f{j}).est.P;
            rec.agent{i}.full.(f{j}).xhat{m}(:,k) = agent{i}.full.(f{j}).est.xhat;
            rec.agent{i}.full.(f{j}).P(:,:,k) = agent{i}.full.(f{j}).est.P;
        end
        % REF:
        rec.agent{i}.ref.dca.xhat{m}(:,k) = agent{i}.ref.dca.est.xhat;
        rec.agent{i}.ref.dca.P(:,:,k) = agent{i}.ref.dca.est.P;
        rec.agent{i}.ref.full.xhat{m}(:,k) = agent{i}.ref.full.est.xhat;
        rec.agent{i}.ref.full.P(:,:,k) = agent{i}.ref.full.est.P;
    end
    rec.xt{m} = xt;
end



% --- COMMUNICATION REDUCTION ---------------------------------------------
function Psi = gevokf_dr(R1,R2,m)
    Q = R1*R1;
    S = R1+R2;
    [X,G] = eig(Q,S);
    idx = get_max_idx_vec(G,m);
    Psi = gram_schmidt_procedure(X(:,idx)');
end

function Psi = gevoci_dr(R1,R2,m)
    max_iter = 10; epsilon = 1e-4;
    wmin = 0.001; wmax = 1;
    I1 = inv(R1); 
    w = 0.5; J = Inf;
    
    for k = 1:max_iter
    
        J_prev = J;
    
        A = R1/w; B = R2/(1-w);
        Q = A^2; S = A+B;
        [X,G] = eig(Q,S);
        idx = get_max_idx_vec(G,m);
        Psi = X(:,idx)';
        RPsi = Psi*R2*Psi';
        
        f = @(w) trace( inv(w*I1 + (1-w)*Psi'/RPsi*Psi) );
        w = fminbnd(f,wmin,wmax);
        J = gevoci_loss_function(I1,RPsi,Psi,w);
    
        if abs(J_prev-J)/J < epsilon; break; end
    end
end

function Psi = gevole_dr(R1,R2,m)
    nx = size(R1,1); 
    [U1,D1] = eig(make_symmetric(inv(R1))); 
    T1 = inv(sqrtm(D1))*U1';
    [U2,D2] = eig(make_symmetric(T1/R2*T1'));
    T = U2'*T1; Ti = inv(T);
    Gammai = Ti*diag(min([ones(1,nx) ; diag(D2)']))*Ti';
    R12 = R1*Gammai*R2; R12 = 0.99*R12;
    Q = make_symmetric((R1-R12)'*(R1-R12)); 
    S = make_symmetric(R1+R2-R12-R12');
    [X,G] = eig(Q,S,'qz');
    idx = get_max_idx_vec(G,m);
    Psi = gram_schmidt_process(X(:,idx)');
end

function J = gevoci_loss_function(I1,RPsi,Psi,w)
    J = trace( inv(w*I1+(1-w)*Psi'/RPsi*Psi) );
end

function U = gram_schmidt_procedure(X)
    [n,m] = size(X);
    is_tall = n > m;    
    if ~is_tall; X = X'; [n,m] = size(X); end 
    
    U = zeros(n,m);
    U(:,1) = X(:,1)/norm(X(:,1));
    
    for i = 2:m
        x = X(:,i); u = x;
        for j = 1:i-1
            uprev = U(:,j);
            u = u - po(uprev,x);
        end
        U(:,i) = u/norm(u);
    end    
    if ~is_tall; U = U'; end  
end

function vu = po(u,v)
    vu = (u'*v/(u'*u))*u;
end

function DS2 = dca_eig(R2)
    D2 = diag(diag(R2));
    T = sqrtm(inv(D2));
    s = max(eig(T*R2*T));
    DS2 = s*D2;
end



% --- FUSION METHODS ------------------------------------------------------
function [xhat,P] = kf_fusion(y1,R1,y2,R2,H2)
    I1 = inv(R1);
    P = inv(I1 + H2'/R2*H2);
    xhat = P*(I1*y1 + H2'/R2*y2);
end

function [xhat,P] = ci_fusion(y1,R1,y2,R2,H2)
    I1 = inv(R1); 
    J = @(w) trace(inv(w*I1 + (1-w)*H2'/R2*H2));
    w = fminbnd(J,0,1);
    P = inv(w*I1 + (1-w)*H2'/R2*H2);
    xhat = P*(w*I1*y1 + (1-w)*H2'/R2*y2);
end

function [xhat,P] = le_fusion(y1,R1,y2,R2,H2)
    nx = length(y1);
    I1 = inv(R1); 
    I2 = make_symmetric(H2'/R2*H2);
    [U1,D1] = eig(make_symmetric(I1));
    T1 = inv(sqrtm(D1))*U1.';
    [U2,D2] = eig(make_symmetric(T1*I2*T1'));
    T = U2.'*T1; Ti = inv(T);
    iota1 = T*I1*y1; iota2 = T*H2'/R2*y2;

    iota = zeros(nx,1);
    D = zeros(nx);
    for i = 1:nx
        if D2(i,i) < 1; iota(i) = iota1(i); D(i,i) = 1;
        else; iota(i) = iota2(i); D(i,i) = D2(i,i);
        end
    end
    P = make_symmetric(inv(Ti*D*Ti'));
    xhat = P*Ti*iota;
end




% --- BASIC FUNCTIONS -----------------------------------------------------
function xt = simulate_target_trajectory(par,process)
    nx = size(process.F,1);
    xt = zeros(nx,par.nk);
    xt(:,1) = mvnrnd(zeros(1,nx),process.P0)';
    for k = 2:par.nk
        w = mvnrnd(zeros(1,nx),process.Q)';
        xt(:,k) = process.F*xt(:,k-1) + w;
    end
end

function A = get_agent_struct(par)
    c = cell(1,par.nk);
    s.K = c; s.P = c; s.Ptu = c; s.Pmu = c; s.est.xhat = []; s.est.P = []; s.est.H = [];
    f = get_method_fieldnames();
    for i = 1:length(f) 
        A.dr.(f{i}) = s; 
        A.full.(f{i}) = s; 
    end
    A.ref.dca = s;
    A.ref.full = s;
    A.sensor.H = [eye(2) zeros(2)]; A.sensor.C = NaN;
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
    T = zeros(par.nagents);
    T(par.nagents,1) = 1; 
    for i = 1:par.nagents-1; T(i,i+1) = 1; end
end

function fields = get_method_fieldnames()
    fields = {'kf' ; 'ci' ; 'le'};
end

function rec = get_datarecord(par,process)
    nx = size(process.F,1);
    c = cell(par.M,1); a = zeros(nx,par.nk);
    for m = 1:par.M; c{m} = a; end
    rec.agent = cell(par.nagents,1);
    rec.xt = c;
    f = get_method_fieldnames();
    for i = 1:length(f); s.(f{i}).xhat = c; end
    s.ref.dca.xhat = c; 
    for j = 1:par.nagents; rec.agent{j} = s; end
end



