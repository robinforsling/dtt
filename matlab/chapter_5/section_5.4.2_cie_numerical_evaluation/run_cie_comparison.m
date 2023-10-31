function data = run_cie_comparison(par,m_vec)
% --- run_cie_comparison() ------------------------------------------------
% Comparison between using CIE and global knowledge.
%
% 2023-10-30 Robin Forsling

p_idx = 1:2; 
SEED = par.SEED;


% --- INITIALIZE ---
scen = load_scen(par); 
nagents = scen.nagents; ntargets = scen.ntargets;

itgt = 1;
nm = length(m_vec);

% GEVO
data.rkf_loc = cell(nm,nagents); data.pkf_loc = data.rkf_loc; data.akf_loc = data.rkf_loc;  
data.rkf_glob = cell(nm,nagents); data.pkf_glob = data.rkf_glob; data.akf_glob = data.rkf_glob;  

data.rci_loc = cell(nm,nagents); data.pci_loc = data.rci_loc; data.aci_loc = data.rci_loc; 
data.rci_glob = cell(nm,nagents); data.pci_glob = data.rci_glob; data.aci_glob = data.rci_glob;

data.rle_loc = cell(nm,nagents); data.ple_loc = data.rle_loc; data.ale_loc = data.rle_loc; 
data.rle_glob = cell(nm,nagents); data.ple_glob = data.rle_glob; data.ale_glob = data.rle_glob;

% DCA
data.rdeci = cell(1,nagents); data.pdeci = data.rdeci; data.adeci = data.rdeci;

data.rlkf = cell(1,nagents); data.plkf = data.rlkf; data.alkf = data.rlkf; 



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
        [data.rkf_loc{im,ia},data.pkf_loc{im,ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia},p_idx);
        [~,data.akf_loc{im,ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end
    
    % --- CI ---
    fprintf('\nRunning CI...\n')
    par.comm_mgmt_method = def_communication_management_technique('gevo');
    par.dimred_loss = def_dimred_loss_function('ci');
    par.fus_method = def_fusion_method('ci');
    rng(SEED) 
    sim_out = simenv(par);  
    for ia = 1:nagents
        [data.rci_loc{im,ia},data.pci_loc{im,ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia},p_idx);
        [~,data.aci_loc{im,ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end

    % --- LE ---
    fprintf('\nRunning LE...\n')
    par.comm_mgmt_method = def_communication_management_technique('gevo');
    par.dimred_loss = def_dimred_loss_function('le');
    par.fus_method = def_fusion_method('le');
    rng(SEED) 
    sim_out = simenv(par);  
    for ia = 1:nagents
        [data.rle_loc{im,ia},data.ple_loc{im,ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia},p_idx);
        [~,data.ale_loc{im,ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
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
        [data.rkf_glob{im,ia},data.pkf_glob{im,ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia},p_idx);
        [~,data.akf_glob{im,ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end
    
    % --- CI ---
    fprintf('\nRunning CI...\n')
    par.comm_mgmt_method = def_communication_management_technique('gevo');
    par.dimred_loss = def_dimred_loss_function('ci');
    par.fus_method = def_fusion_method('ci');
    rng(SEED) 
    sim_out = simenv(par); 
    for ia = 1:nagents
        [data.rci_glob{im,ia},data.pci_glob{im,ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia},p_idx);
        [~,data.aci_glob{im,ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end

    % --- LE ---
    fprintf('\nRunning LE...\n')
    par.comm_mgmt_method = def_communication_management_technique('gevo');
    par.dimred_loss = def_dimred_loss_function('le');
    par.fus_method = def_fusion_method('le');
    rng(SEED) 
    sim_out = simenv(par);  
    for ia = 1:nagents
        [data.rle_glob{im,ia},data.ple_glob{im,ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia},p_idx);
        [~,data.ale_glob{im,ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end
end

fprintf('\n--- DCA ---\n')
fprintf('\nRunning DCA-EIG...\n')
par.comm_mgmt_method = def_communication_management_technique('dca-eig');
par.fus_method = def_fusion_method('ci');
rng(SEED) 
sim_out = simenv(par);
for ia = 1:nagents
    [data.rdeci{ia},data.pdeci{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
    [~,data.adeci{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
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
    [data.rlkf{ia},data.plkf{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
    [~,data.alkf{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
end

% --- CRLB ---
fprintf('\nComputing CRLB...\n')
nagents = sim_out.scen.nagents;
sm = cell(nagents,1);
for ia = 1:nagents
    sm{ia} = sim_out.scen.agents{ia}.sensor_model;
end
pm = sim_out.scen.agents{1}.process_model; 
[data.PCRLB,data.pcrlb] = crlb(sim_out.scen.XT{itgt},sim_out.scen.XS,sm,pm); 

t_elapsed = toc(t_start); t_end_str = get_datetime;
fprintf('\n\nSimulations started at %s\n',t_start_str)
fprintf('Simulations finished at %s\n',t_end_str)
fprintf('Elapsed time: %s\n',get_elapsed_time(t_elapsed))

data.nx = size(sim_out.rec{1}.xhat{1},1);
data.N = scen.par.N;
data.M = length(sim_out.rec{1}.P);

end

function str = get_elapsed_time(te)
    hours = floor(te/3600);
    te = te - 3600*hours;
    minutes = floor(te/60);
    seconds = te - 60*minutes;
    str = sprintf('%d:%2.0f:%2.0f',hours,minutes,round(seconds));
end
