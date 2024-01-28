function [out_par,data,crlb] = run_cie_comparison(par,m_vec)
% --- run_cie_comparison() ------------------------------------------------
% Comparison between using CIE and global knowledge.
%
% 2023-10-30 Robin Forsling

p_idx = 1:2; 
SEED = par.SEED;


% --- INITIALIZE ---
scen = load_scen(par); 
nagents = scen.nagents; 
itgt = 1;
nm = length(m_vec);

data = cell(1,nagents);
for ia = 1:nagents
    data{ia}.kf = cell(nm,1);
    data{ia}.ci = cell(nm,1);
    data{ia}.le = cell(nm,1);
end



% --- SIMULATIONS ---
t_start = tic; t_start_str = get_datetime;
fprintf('\n\nSimulations starting at %s\n',t_start_str)
fprintf('\n--- GEVO: LOCAL INFORMATION ---\n')

par.cntrl.globally_known_est = 0;
for im = 1:nm

    par.m = m_vec(im);
    fprintf('\n--- m = %d ---\n',par.m)

    % --- KF ---
    fprintf('\nRunning KF...\n')
    par.comm_mgmt_method = def_communication_management_technique('gevo');
    par.dimred_loss = def_dimred_loss_function('kf-rce');
    par.fus_method = def_fusion_method('gimf');
    rng(SEED)
    sim_out = simenv(par);    
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        data{ia}.kf{im}.loc = merge_structs(c,p);
    end
    
    % --- CI ---
    fprintf('\nRunning CI...\n')
    par.comm_mgmt_method = def_communication_management_technique('gevo');
    par.dimred_loss = def_dimred_loss_function('ci');
    par.fus_method = def_fusion_method('ci');
    rng(SEED) 
    sim_out = simenv(par);  
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        data{ia}.ci{im}.loc = merge_structs(c,p);
    end

    % --- LE ---
    fprintf('\nRunning LE...\n')
    par.comm_mgmt_method = def_communication_management_technique('gevo');
    par.dimred_loss = def_dimred_loss_function('le');
    par.fus_method = def_fusion_method('le');
    rng(SEED) 
    sim_out = simenv(par);  
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        data{ia}.le{im}.loc = merge_structs(c,p);
    end
end

fprintf('\n--- GEVO: GLOBAL INFORMATION ---\n')

par.cntrl.globally_known_est = 1;
for im = 1:nm

    par.m = m_vec(im);
    fprintf('\n--- m = %d ---\n',par.m)

    % --- KF ---
    fprintf('\nRunning KF...\n')
    par.comm_mgmt_method = def_communication_management_technique('gevo');
    par.dimred_loss = def_dimred_loss_function('kf-rce');
    par.fus_method = def_fusion_method('gimf');
    rng(SEED)
    sim_out = simenv(par);
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        data{ia}.kf{im}.glob = merge_structs(c,p);
    end
    
    % --- CI ---
    fprintf('\nRunning CI...\n')
    par.comm_mgmt_method = def_communication_management_technique('gevo');
    par.dimred_loss = def_dimred_loss_function('ci');
    par.fus_method = def_fusion_method('ci');
    rng(SEED) 
    sim_out = simenv(par); 
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        data{ia}.ci{im}.glob = merge_structs(c,p);
    end

    % --- LE ---
    fprintf('\nRunning LE...\n')
    par.comm_mgmt_method = def_communication_management_technique('gevo');
    par.dimred_loss = def_dimred_loss_function('le');
    par.fus_method = def_fusion_method('le');
    rng(SEED) 
    sim_out = simenv(par);  
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        data{ia}.le{im}.glob = merge_structs(c,p);
    end
end

fprintf('\n--- DCA ---\n')
fprintf('\nRunning DCA-EIG...\n')
par.comm_mgmt_method = def_communication_management_technique('dca-eig');
par.fus_method = def_fusion_method('ci');
rng(SEED) 
sim_out = simenv(par);
for ia = 1:nagents
    X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
    c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
    data{ia}.dca = merge_structs(c,p);
end

fprintf('\n--- MISC ---\n')

% --- LKF ---
fprintf('\nRunning LKF...\n')    
par.comm_mgmt_method = def_communication_management_technique('no-communication');
par.dimred_loss = 0;
par.fus_method = def_fusion_method('lkf');
par.cntrl.globally_known_est = 0;  
rng(SEED)
sim_out = simenv(par);
for ia = 1:nagents
    X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
    c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
    data{ia}.lkf = merge_structs(c,p);
end

% --- CRLB ---
fprintf('\nComputing CRLB...\n')
nagents = sim_out.nagents;
sm = cell(nagents,1);
for ia = 1:nagents
    sm{ia} = sim_out.agents{ia}.sensor_model;
end
pm = sim_out.agents{1}.process_model;   
X = sim_out.targets{itgt}.X; 
crlb = cramer_rao_lower_bound(X,sim_out.sensor_pos,sm,pm); 

t_elapsed = toc(t_start); t_end_str = get_datetime;
fprintf('\n\nSimulations started at %s\n',t_start_str)
fprintf('Simulations finished at %s\n',t_end_str)
fprintf('Elapsed time: %s\n',get_elapsed_time(t_elapsed))

out_par.nx = size(sim_out.rec{1}.xhat{1},1);
out_par.N = scen.par.N;
out_par.M = length(sim_out.rec{1}.P);

end

function str = get_elapsed_time(te)
    hours = floor(te/3600);
    te = te - 3600*hours;
    minutes = floor(te/60);
    seconds = te - 60*minutes;
    str = sprintf('%d:%2.0f:%2.0f',hours,minutes,round(seconds));
end
