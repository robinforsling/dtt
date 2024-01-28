function dca_evaluation
% --- dca_evaluation() ----------------------------------------------------
% Section 3.3.4 Numerical Evaluation - DCA: Figure 3.11
%
% 2024-01-16 Robin Forsling

set_latex_interpreter;

SEED = round(1000*rand);


% TOGGLE
KF_ON = 1; 
CI_ON = 1; 
LKF_ON = 1;
NKF_ON = 1;
FULL_ON = 1;
CRLB_ON = 1;


% PARAMETERS
par = default_simulation_parameters;
par.ID = 1;
par.sigma_r = 1000;
par.sigma_az = 1.0*d2r;
par.M = 10;


% --- SIMULATIONS ---
scen = load_scen(par); 
nagents = scen.nagents;
itgt = 1; 
eval_res = cell(1,nagents);

fprintf('\n\nSimulations starting at %s\n',get_datetime)
fprintf('\n--- DCA ---\n')

% --- CI ---
if CI_ON
    fprintf('\nRunning CI...\n')

    par.dimred_loss = def_dimred_loss_function('ci');

    % DCA-EIG
    fprintf('\nDCA-EIG:\n')
    par.comm_mgmt_method = def_communication_management_technique('dca-eig');
    par.fus_method = def_fusion_method('ci');
    rng(SEED) 
    sim_out = simenv(par);
    
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        eval_res{ia}.deci = merge_structs(c,p);
    end

    % DCA-OPT
    fprintf('\nDCA-OPT:\n')
    par.comm_mgmt_method = def_communication_management_technique('dca-opt');
    par.fus_method = def_fusion_method('ci');
    rng(SEED) 
    sim_out = simenv(par);
    
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        eval_res{ia}.doci = merge_structs(c,p);
    end

    % DCA-DOM
    fprintf('\nDCA-DOM:\n')
    par.comm_mgmt_method = def_communication_management_technique('dca-dom');
    par.fus_method = def_fusion_method('ci');
    rng(SEED) 
    sim_out = simenv(par);
    
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        eval_res{ia}.ddci = merge_structs(c,p);
    end

    % DCA-HYP
    fprintf('\nDCA-HYP:\n')
    par.comm_mgmt_method = def_communication_management_technique('dca-hyp');
    par.fus_method = def_fusion_method('ci-dca-hyp');
    rng(SEED) 
    sim_out = simenv(par); 
    
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        eval_res{ia}.dhci = merge_structs(c,p);
    end
end
    


if FULL_ON
    fprintf('\n--- FULL ---\n')

    if CI_ON 
        fprintf('\nRunning CI...\n')
        par.comm_mgmt_method = def_communication_management_technique('full');
        par.fus_method = def_fusion_method('ci');
        rng(SEED)
        sim_out = simenv(par);

        for ia = 1:nagents
            X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
            c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
            eval_res{ia}.fci = merge_structs(c,p);
        end
    end
end

fprintf('\n--- MISC ---\n')

% --- NKF ---
if NKF_ON
    fprintf('\nRunning NKF...\n')    
    par.comm_mgmt_method = def_communication_management_technique('full');
    par.fus_method = def_fusion_method('kf');
    rng(SEED)
    sim_out = simenv(par);
    
    for ia = 1:nagents
        X = sim_out.targets{itgt}.X; rec = sim_out.rec{ia};
        c = credibility_assessment(X,rec); p = performance_assessment(X,rec);
        eval_res{ia}.nkf = merge_structs(c,p);
    end
end

% --- LKF ---
if LKF_ON
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
        eval_res{ia}.lkf = merge_structs(c,p);
    end
end

% --- CRLB ---
if CRLB_ON
    fprintf('\nComputing CRLB...\n')
    nagents = sim_out.nagents;
    sm = cell(nagents,1);
    for i = 1:nagents; sm{i} = sim_out.agents{i}.sensor_model; end
    pm = sim_out.agents{1}.process_model; pm.q = pm.q; % :D    
    X = sim_out.targets{itgt}.X; 
    crlb = cramer_rao_lower_bound(X,sim_out.sensor_pos,sm,pm); 
end

fprintf('\nSimulations finished at %s\n',get_datetime)




%% --- PLOT ---
nrow = 3;
ncol = nagents;

c = 1/par.sigma_r;

n = size(sim_out.rec{1}.xhat{1},1);
M = length(sim_out.rec{1}.P);
anees_ci = anees_confidence_interval(n,M);

N = scen.par.N; t = 1:N;
idx = 1:1:N;

lw = 1.5; lwb = 2.0; lwf = 2.0;
clr = get_thesis_colors;
clrde = clr.orange; clrdo = clr.darkyellow; clrdd = clr.green;
clrdh = clr.blue; clrf = clr.purple; clr.lkf = clr.darkgray;


figure(1);clf
for ia = 1:nagents

    res = eval_res{ia};


    % --- RMT ---
    subplot(nrow,ncol,ia); hold on

    if CRLB_ON; h = plot(t(2:end),c*crlb.glob.rt.pos(2:end),'k-','DisplayName','CRLB'); h.LineWidth = lwb; end
    if LKF_ON; h = plot(t(idx),c*res.lkf.rmt.pos(idx),'k-','DisplayName','LKF'); h.Color = clr.lkf; h.LineWidth = lwb; end
    if NKF_ON; h = plot(t(idx),c*res.nkf.rmt.pos(idx),'k:','DisplayName','NKF'); h.Color = clr.naive; h.LineWidth = lwb; end
    
    if CI_ON
        h = plot(t(idx),c*res.deci.rmt.pos(idx),'-','DisplayName','DCA-EIG-CI'); h.LineWidth = lw; h.Color = clrde;
        h = plot(t(idx),c*res.doci.rmt.pos(idx),'-','DisplayName','DCA-OPT-CI'); h.LineWidth = lw; h.Color = clrdo;
        h = plot(t(idx),c*res.ddci.rmt.pos(idx),'-','DisplayName','DCA-DOM-CI'); h.LineWidth = lw; h.Color = clrdd;
        h = plot(t(idx),c*res.dhci.rmt.pos(idx),'-','DisplayName','DCA-HYP-CI'); h.LineWidth = lw; h.Color = clrdh;
        if FULL_ON; h = plot(t(idx),c*res.fci.rmt.pos(idx),'--','DisplayName','FULL-CI'); h.LineWidth = lwf; h.Color = clrf; end
    end

    ylim([0 c*res.lkf.rmt.pos(1)])
    xlabel('$k$','interpreter','latex'); 
    ylabel('RMT$/\sigma_r$','interpreter','latex')
    title(sprintf('Agent %d',ia))
    if ia == 1; legend show; end
    box on


    % --- RMSE ---
    subplot(nrow,ncol,ia+ncol); hold on

    if CRLB_ON; h = plot(t(2:end),c*crlb.glob.rt.pos(2:end),'k-','DisplayName','CRLB'); h.LineWidth = lwb; end
    if LKF_ON; h = plot(t(idx),c*res.lkf.rmse.pos(idx),'k-','DisplayName','LKF'); h.Color = clr.lkf; h.LineWidth = lwb; end
    if NKF_ON; h = plot(t(idx),c*res.nkf.rmse.pos(idx),'k:','DisplayName','NKF'); h.Color = clr.naive; h.LineWidth = lwb; end
    
    if CI_ON
        h = plot(t(idx),c*res.deci.rmse.pos(idx),'-','DisplayName','DCA-EIG-CI'); h.LineWidth = lw; h.Color = clrde;
        h = plot(t(idx),c*res.doci.rmse.pos(idx),'-','DisplayName','DCA-OPT-CI'); h.LineWidth = lw; h.Color = clrdo;
        h = plot(t(idx),c*res.ddci.rmse.pos(idx),'-','DisplayName','DCA-DOM-CI'); h.LineWidth = lw; h.Color = clrdd;
        h = plot(t(idx),c*res.dhci.rmse.pos(idx),'-','DisplayName','DCA-HYP-CI'); h.LineWidth = lw; h.Color = clrdh;
        if FULL_ON; h = plot(t(idx),c*res.fci.rmse.pos(idx),'--','DisplayName','FULL-CI'); h.LineWidth = lwf; h.Color = clrf; end
    end

    ylim([0 c*res.lkf.rmse.pos(1)])
    xlabel('$k$','interpreter','latex'); 
    ylabel('RMSE$/\sigma_r$','interpreter','latex')
    box on


    % --- ANEES ---
    x = t(idx); x2 = [x fliplr(x)];
    subplot(nrow,ncol,ia+2*ncol); hold on

    h = plot(t(idx),anees_ci.lvl999(1)*ones(1,length(idx)),'k--'); turn_legend_item_off(h);
    h = plot(t(idx),anees_ci.lvl999(2)*ones(1,length(idx)),'k--','DisplayName','99\% conf. int.'); 

    if LKF_ON; h = plot(t(idx),res.lkf.anees(idx),'k-','DisplayName','LKF'); h.Color = clr.lkf; h.LineWidth = lwb; end
    if NKF_ON; h = plot(t(idx),res.nkf.anees(idx),'k:','DisplayName','NKF'); h.Color = clr.naive; h.LineWidth = lwb; end

    if CI_ON
        h = plot(t(idx),res.deci.anees(idx),'-','DisplayName','DCA-EIG-CI'); h.LineWidth = lw; h.Color = clrde; 
        h = plot(t(idx),res.doci.anees(idx),'-','DisplayName','DCA-OPT-CI'); h.LineWidth = lw; h.Color = clrdo;
        h = plot(t(idx),res.ddci.anees(idx),'-','DisplayName','DCA-DOM-CI'); h.LineWidth = lw; h.Color = clrdd;
        h = plot(t(idx),res.dhci.anees(idx),'-','DisplayName','DCA-HYP-CI'); h.LineWidth = lw; h.Color = clrdh;
        if FULL_ON; h = plot(t(idx),res.fci.anees(idx),'--','DisplayName','FULL-CI'); h.LineWidth = lwf; h.Color = clrf; end
    end

    ylim([0.25 2.25])
    xlabel('$k$','interpreter','latex'); ylabel('ANEES','interpreter','latex')
    if ia == 1; legend show; end
    box on

    set_fontsize_all(14)

end

