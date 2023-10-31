function dca_evaluation
% --- dca_evaluation() ----------------------------------------------------
% Section 3.3.4 Numerical Evaluation - DCA: Figure 3.11
%
% 2023-10-30 Robin Forsling

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
par = get_default_sim_params;
par.ID = 1;
par.sigma_r = 1000;
par.sigma_az = 1.0*d2r;
par.M = 100;


% --- SIMULATIONS ---
scen = load_scen(par); 
nagents = scen.nagents; ntargets = scen.ntargets;

itgt = 1; itgt = min([itgt ntargets]);

if CI_ON 
    rdeci = cell(1,nagents); pdeci = rdeci; adeci = rdeci; 
    rdoci = cell(1,nagents); pdoci = rdoci; adoci = rdeci;
    rddci = cell(1,nagents); pddci = rddci; addci = rddci;
    rdhci = cell(1,nagents); pdhci = rdhci; adhci = rdhci;
end

if FULL_ON
    if KF_ON; rfkf = cell(1,nagents); pfkf = rfkf; afkf = rfkf; end
    if CI_ON; rfci = cell(1,nagents); pfci = rfci; afci = rfci; end
end

if LKF_ON; rlkf = cell(1,nagents); plkf = rlkf; alkf = rlkf; end
if NKF_ON; rnkf = cell(1,nagents); pnkf = rnkf; ankf = rnkf; end

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
        [rdeci{ia},pdeci{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
        [~,adeci{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end

    % DCA-OPT
    fprintf('\nDCA-OPT:\n')
    par.comm_mgmt_method = def_communication_management_technique('dca-opt');
    par.fus_method = def_fusion_method('ci');
    rng(SEED) 
    sim_out = simenv(par);
    
    for ia = 1:nagents
        [rdoci{ia},pdoci{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
        [~,adoci{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end

    % DCA-DOM
    fprintf('\nDCA-DOM:\n')
    par.comm_mgmt_method = def_communication_management_technique('dca-dom');
    par.fus_method = def_fusion_method('ci');
    rng(SEED) 
    sim_out = simenv(par);
    
    for ia = 1:nagents
        [rddci{ia},pddci{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
        [~,addci{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end

    % DCA-HYP
    fprintf('\nDCA-HYP:\n')
    par.comm_mgmt_method = def_communication_management_technique('dca-hyp');
    par.fus_method = def_fusion_method('ci-dca-hyp');
    rng(SEED) 
    sim_out = simenv(par); 
    
    for ia = 1:nagents
        [rdhci{ia},pdhci{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
        [~,adhci{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
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
            [rfci{ia},pfci{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
            [~,afci{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
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
        [rnkf{ia},pnkf{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
        [~,ankf{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
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
        [rlkf{ia},plkf{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
        [~,alkf{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end
end

% --- CRLB ---
if CRLB_ON
    fprintf('\nComputing CRLB...\n')
    nagents = sim_out.scen.nagents;
    sm = cell(nagents,1);
    for i = 1:nagents
        sm{i} = sim_out.scen.agents{i}.sensor_model;
    end
    pm = sim_out.scen.agents{1}.process_model;
    pm.q = pm.q; % :D
    
    [PCRLB,pcrlb] = crlb(sim_out.scen.XT{itgt},sim_out.scen.XS,sm,pm); 
end

fprintf('\nSimulations finished at %s\n',get_datetime)




% --- PLOT ---
nrow = 3;
ncol = nagents;

c = 1/par.sigma_r;

n = size(sim_out.rec{1}.xhat{1},1);
M = length(sim_out.rec{1}.P);
anees_ci = anees_confidence_interval(n,M);

N = scen.par.N; t = 1:N;
idx = 1:1:N;
lw = 1.5;
lwb = 2.0;
lwf = 2.0;
clr = get_thesis_colors;
clrde = clr.orange;
clrdo = clr.darkyellow; 
clrdd = clr.green;
clrdh = clr.blue;
clrf = clr.purple;
clr.lkf = clr.darkgray;


figure(1);clf
for ia = 1:nagents


    % --- RMT ---
    subplot(nrow,ncol,ia); hold on

    if CRLB_ON; h = plot(t(2:end),c*pcrlb(2:end),'k-','DisplayName','CRLB'); h.LineWidth = lwb; end
    if LKF_ON; h = plot(t(idx),c*plkf{ia}(idx),'k-','DisplayName','LKF'); h.Color = clr.lkf; h.LineWidth = lwb; end
    if NKF_ON; h = plot(t(idx),c*pnkf{ia}(idx),'k:','DisplayName','NKF'); h.Color = clr.naive; h.LineWidth = lwb; end
    
    if CI_ON
        h = plot(t(idx),c*pdeci{ia}(idx),'-','DisplayName','DCA-EIG-CI'); h.LineWidth = lw; h.Color = clrde;
        h = plot(t(idx),c*pdoci{ia}(idx),'-','DisplayName','DCA-OPT-CI'); h.LineWidth = lw; h.Color = clrdo;
        h = plot(t(idx),c*pddci{ia}(idx),'-','DisplayName','DCA-DOM-CI'); h.LineWidth = lw; h.Color = clrdd;
        h = plot(t(idx),c*pdhci{ia}(idx),'-','DisplayName','DCA-HYP-CI'); h.LineWidth = lw; h.Color = clrdh;
        if FULL_ON; h = plot(t(idx),c*pfci{ia}(idx),'--','DisplayName','FULL-CI'); h.LineWidth = lwf; h.Color = clrf; end
    end

    ylim([0 c*plkf{ia}(1)])
    xlabel('$k$','interpreter','latex'); 
    ylabel('RMT$/\sigma_r$','interpreter','latex')
    title(sprintf('Agent %d',ia))
    if ia == 1; legend show; end
    box on


    % --- RMSE ---
    subplot(nrow,ncol,ia+ncol); hold on

    if CRLB_ON; h = plot(t(2:end),c*pcrlb(2:end),'k-','DisplayName','CRLB'); h.LineWidth = lwb; end
    if LKF_ON; h = plot(t(idx),c*rlkf{ia}(idx),'k-','DisplayName','LKF'); h.Color = clr.lkf; h.LineWidth = lwb; end
    if NKF_ON; h = plot(t(idx),c*rnkf{ia}(idx),'k:','DisplayName','NKF'); h.Color = clr.naive; h.LineWidth = lwb; end
    
    if CI_ON
        h = plot(t(idx),c*rdeci{ia}(idx),'-','DisplayName','DCA-EIG-CI'); h.LineWidth = lw; h.Color = clrde;
        h = plot(t(idx),c*rdoci{ia}(idx),'-','DisplayName','DCA-OPT-CI'); h.LineWidth = lw; h.Color = clrdo;
        h = plot(t(idx),c*rddci{ia}(idx),'-','DisplayName','DCA-DOM-CI'); h.LineWidth = lw; h.Color = clrdd;
        h = plot(t(idx),c*rdhci{ia}(idx),'-','DisplayName','DCA-HYP-CI'); h.LineWidth = lw; h.Color = clrdh;
        if FULL_ON; h = plot(t(idx),c*rfci{ia}(idx),'--','DisplayName','FULL-CI'); h.LineWidth = lwf; h.Color = clrf; end
    end

    ylim([0 c*rlkf{ia}(1)])
    xlabel('$k$','interpreter','latex'); 
    ylabel('RMSE$/\sigma_r$','interpreter','latex')
    box on


    % --- ANEES ---
    x = t(idx); x2 = [x fliplr(x)];
    subplot(nrow,ncol,ia+2*ncol); hold on

    h = plot(t(idx),anees_ci.lvl999(1)*ones(1,length(idx)),'k--'); turn_legend_item_off(h);
    h = plot(t(idx),anees_ci.lvl999(2)*ones(1,length(idx)),'k--','DisplayName','99\% conf. int.'); 

    if LKF_ON; h = plot(t(idx),alkf{ia}(idx),'k-','DisplayName','LKF'); h.Color = clr.lkf; h.LineWidth = lwb; end
    if NKF_ON; h = plot(t(idx),ankf{ia}(idx),'k:','DisplayName','NKF'); h.Color = clr.naive; h.LineWidth = lwb; end

    if CI_ON
        h = plot(t(idx),adeci{ia}(idx),'-','DisplayName','DCA-EIG-CI'); h.LineWidth = lw; h.Color = clrde; 
        h = plot(t(idx),adoci{ia}(idx),'-','DisplayName','DCA-OPT-CI'); h.LineWidth = lw; h.Color = clrdo;
        h = plot(t(idx),addci{ia}(idx),'-','DisplayName','DCA-DOM-CI'); h.LineWidth = lw; h.Color = clrdd;
        h = plot(t(idx),adhci{ia}(idx),'-','DisplayName','DCA-HYP-CI'); h.LineWidth = lw; h.Color = clrdh;
        if FULL_ON; h = plot(t(idx),afci{ia}(idx),'--','DisplayName','FULL-CI'); h.LineWidth = lwf; h.Color = clrf; end
    end

    ylim([0.25 2.25])
    xlabel('$k$','interpreter','latex'); ylabel('ANEES','interpreter','latex')
    if ia == 1; legend show; end
    box on

    set_fontsize_all(14)

end

