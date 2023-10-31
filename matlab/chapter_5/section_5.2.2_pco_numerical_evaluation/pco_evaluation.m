function pco_evaluation
% --- pco_evaluation() ----------------------------------------------------
% Secion 5.2.2 Numerical Evaluation - PCO: Figure 5.4
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- PARAMETERS ----------------------------------------------------------
par.m = 2;
par.nagents = 3;
par.nk = 15;
par.M = 100;
par.network_connectivity = 1;                                               % 1 = reduced, 2 = full
                                                        

% --- AGENTS --------------------------------------------------------------
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


% --- ESTIMATION LOOP -----------------------------------------------------
fprintf('\n--- PCO EVALUATION ---\n\n')
fprintf('m = %d\n',par.m)
fprintf('nagents = %d\n',par.nagents)
fprintf('nk = %d\n',par.nk)
fprintf('M = %d\n\n',par.M)

switch par.network_connectivity
    case 1; fprintf('Network connnectivity: reduced\n')
    case 2; fprintf('Network connnectivity: full\n')
end

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

        % PCO & COMMUNICATION:
        data_tx = manage_communication(par,agent{idx_tx{k}});

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


% --- GENERATE TIKZ CODE --------------------------------------------------
%generate_tikz_code(rec,stats,par)

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
                P = rec_ag.P(:,:,k);
                S.(f{i}){k} = E;
                coin.(f{i})(k) = compute_coin(P,E);
                anees.(f{i})(k) = trace(E/P)/nx;
                rmt.(f{i})(k) = sqrt(trace(P));
                rmse.(f{i})(k) = sqrt(trace(E));
            end
        end
        stats{iagent}.dr.coin = coin;
        stats{iagent}.dr.anees = anees;
        stats{iagent}.dr.rmt = rmt;
        stats{iagent}.dr.rmse = rmse;
        stats{iagent}.dr.S = S;
    end

    % FULL:
    for iagent = 1:nagents
        for i = 1:length(f)
            rec_ag = rec.agent{iagent}.full.(f{i});
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
                P = rec_ag.P(:,:,k);
                S.(f{i}){k} = E;
                coin.(f{i})(k) = compute_coin(P,E);
                anees.(f{i})(k) = trace(E/P)/nx;
                rmt.(f{i})(k) = sqrt(trace(P));
                rmse.(f{i})(k) = sqrt(trace(E));
            end
        end
        stats{iagent}.full.coin = coin;
        stats{iagent}.full.anees = anees;
        stats{iagent}.full.rmt = rmt;
        stats{iagent}.full.rmse = rmse;
        stats{iagent}.full.S = S;
    end

    % RATIO:
    for iagent = 1:nagents
        for i = 1:length(f)
            stats{iagent}.ratio.rmt.(f{i}) = stats{iagent}.dr.rmt.(f{i}) ./ stats{iagent}.full.rmt.(f{i});
            stats{iagent}.ratio.rmse.(f{i}) = stats{iagent}.dr.rmse.(f{i}) ./ stats{iagent}.full.rmse.(f{i});
        end
    end
end

function coin = compute_coin(P,E)
    L = chol(P,'lower'); 
    Li = inv(L);
    coin = max(eig(make_symmetric(Li*E*Li')));
end

function plot_stats(stats,iagent,par)
    coin = stats{iagent}.dr.coin;    
    anees = stats{iagent}.dr.anees;
    rmt = stats{iagent}.dr.rmt;
    rmse = stats{iagent}.dr.rmse;
    rmtr = stats{iagent}.ratio.rmt;
    rmser = stats{iagent}.ratio.rmse;

    kvec = 1:par.nk;
    k_fu = get_time_of_fusion(par);
    idx = k_fu{iagent};

    clr = get_thesis_colors;
    lw = 2;

    figure(iagent);clf

    subplot(3,2,1);hold on
    hkf = plot(kvec(idx),coin.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),coin.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hici = plot(kvec(idx),coin.ici(idx),'-'); hici.Color = clr.ici; hici.LineWidth = lw;
    hle = plot(kvec(idx),coin.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('COIN'); box on

    subplot(3,2,2);hold on
    hkf = plot(kvec(idx),anees.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),anees.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hici = plot(kvec(idx),anees.ici(idx),'-'); hici.Color = clr.ici; hici.LineWidth = lw;
    hle = plot(kvec(idx),anees.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('ANEES'); box on

    subplot(3,2,3);hold on
    hkf = plot(kvec(idx),rmt.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),rmt.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hici = plot(kvec(idx),rmt.ici(idx),'-'); hici.Color = clr.ici; hici.LineWidth = lw;
    hle = plot(kvec(idx),rmt.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('RMT'); box on  

    subplot(3,2,4);hold on
    hkf = plot(kvec(idx),rmse.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),rmse.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hici = plot(kvec(idx),rmse.ici(idx),'-'); hici.Color = clr.ici; hici.LineWidth = lw;
    hle = plot(kvec(idx),rmse.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('RMSE'); box on

    subplot(3,2,5);hold on
    hkf = plot(kvec(idx),rmtr.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),rmtr.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hici = plot(kvec(idx),rmtr.ici(idx),'-'); hici.Color = clr.ici; hici.LineWidth = lw;
    hle = plot(kvec(idx),rmtr.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('RMTR'); box on

    subplot(3,2,6);hold on
    hkf = plot(kvec(idx),rmser.kf(idx),'-'); hkf.Color = clr.kf; hkf.LineWidth = lw;
    hci = plot(kvec(idx),rmser.ci(idx),'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hici = plot(kvec(idx),rmser.ici(idx),'-'); hici.Color = clr.ici; hici.LineWidth = lw;
    hle = plot(kvec(idx),rmser.le(idx),'-'); hle.Color = clr.le; hle.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('RMSER'); box on

    legend(gca,[hkf hci hici hle],'KF','CI','ICI','LE')

    set_fontsize_all(14)    

end

function generate_tikz_code(rec,stats,par)

    iagent = 3;
    nt = min(par.M,50);
    f = get_method_fieldnames(); nf = length(f);
    
    kvec = 1:par.nk;
    kfu = get_time_of_fusion(par);
    idx = kfu{iagent};

    coin = stats{iagent}.dr.coin;
    anees = stats{iagent}.dr.anees;
    rmt = stats{iagent}.dr.rmt;
    rmse = stats{iagent}.dr.rmse;
    rmtr = stats{iagent}.ratio.rmt;
    rmser = stats{iagent}.ratio.rmse;

    str_traj = 'data/data__evaluation_pco_trajectory.tex';
    str_data = 'data/data__evaluation_pco_results_';


    % --- TARGET TRAJECTORIES ---
    fileID = fopen(str_traj,'w');
    for m = 1:nt
        str = get_tikz_plot_coordinates(rec.xt{m});
        fprintf(fileID,'\\drawtraj%s\n',str);
    end
    fclose(fileID);


    % --- COIN ---
    fileID = fopen(strcat(str_data,'coin.tex'),'w');
    for i = 1:nf
        str = get_tikz_plot_coordinates(kvec(idx),coin.(f{i})(idx));
        fprintf(fileID,'\\draw%s\\plotcommand%s\n',f{i},str);
    end
    fclose(fileID);


    % --- ANEES ---
    fileID = fopen(strcat(str_data,'anees.tex'),'w');
    for i = 1:nf
        str = get_tikz_plot_coordinates(kvec(idx),anees.(f{i})(idx));
        fprintf(fileID,'\\draw%s\\plotcommand%s\n',f{i},str);
    end
    fclose(fileID);


    % --- RMT ---
    fileID = fopen(strcat(str_data,'rmt.tex'),'w');
    for i = 1:nf
        str = get_tikz_plot_coordinates(kvec(idx),rmt.(f{i})(idx));
        fprintf(fileID,'\\draw%s\\plotcommand%s\n',f{i},str);
    end
    fclose(fileID);


    % --- RMSE ---
    fileID = fopen(strcat(str_data,'rmse.tex'),'w');
    for i = 1:nf
        str = get_tikz_plot_coordinates(kvec(idx),rmse.(f{i})(idx));
        fprintf(fileID,'\\draw%s\\plotcommand%s\n',f{i},str);
    end
    fclose(fileID);


    % --- RMTR ---
    fileID = fopen(strcat(str_data,'rmtr.tex'),'w');
    for i = 1:nf
        str = get_tikz_plot_coordinates(kvec(idx),rmtr.(f{i})(idx));
        fprintf(fileID,'\\draw%s\\plotcommand%s\n',f{i},str);
    end
    fclose(fileID);


    % --- RMSER ---
    fileID = fopen(strcat(str_data,'rmser.tex'),'w');
    for i = 1:nf
        str = get_tikz_plot_coordinates(kvec(idx),rmser.(f{i})(idx));
        fprintf(fileID,'\\draw%s\\plotcommand%s\n',f{i},str);
    end
    fclose(fileID);

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
            case 'ici'; [xhat,P] = ici_fusion(y1,R1,yPsi,RPsi,Psi);
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
            case 'ici'; [xhat,P] = ici_fusion(y1,R1,y2,R2,H2);
            case 'le'; [xhat,P] = le_fusion(y1,R1,y2,R2,H2);
        end
        agent.full.(f{i}).est.xhat = xhat; agent.full.(f{i}).est.P = P;
    end
end

function data_tx = manage_communication(par,agent)
    m = par.m;
    data_tx = agent;
    f = get_method_fieldnames();
    for i = 1:length(f)
        y2 = agent.dr.(f{i}).est.xhat;
        R2 = agent.dr.(f{i}).est.P;
        [U,D] = eig(R2);
        idx = get_min_idx_vec(diag(D),m);
        Psi = U(:,idx)';
        data_tx.dr.(f{i}).est.xhat = Psi*y2;
        data_tx.dr.(f{i}).est.P = Psi*R2*Psi';
        data_tx.dr.(f{i}).est.H = Psi;
    end
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
    end
    rec.xt{m} = xt;
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

function [xhat,P] = ici_fusion(y1,R1,y2,R2,H2)
    I1 = inv(R1);
    J = @(w) trace(inv(I1 + H2'/R2*H2 - w*H2'/(w*H2*R1*H2'+(1-w)*R2)*H2 ...
                    - (1-w)*H2'*H2/(w*R1+(1-w)*H2'*R2*H2)*H2'*H2));
    w = fminbnd(J,0,1);
    G1 = w*H2*R1*H2'+(1-w)*R2; 
    G2 = w*R1+(1-w)*H2'*R2*H2;
    P = inv(I1 + H2'/R2*H2 - w*H2'/G1*H2 - (1-w)*H2'*H2/G2*H2'*H2);
    a = P*(I1 - w*H2'/G1*H2)*y1;
    b = P*(H2'/R2 - (1-w)*H2'*H2/G2*H2')*y2;
    xhat = a + b;
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
        %w = process.LQ*randn(nx,1);
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
    fields = {'kf' ; 'ci' ; 'ici' ; 'le'};
end

function rec = get_datarecord(par,process)
    nx = size(process.F,1);
    c = cell(par.M,1); a = zeros(nx,par.nk);
    for m = 1:par.M; c{m} = a; end
    rec.agent = cell(par.nagents,1);
    rec.xt = c;
    f = get_method_fieldnames();
    for i = 1:length(f); s.(f{i}).xhat = c; end
    for j = 1:par.nagents; rec.agent{j} = s; end
end



