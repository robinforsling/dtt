function run_dtt_scenario

set_latex_interpreter;

SEED = round(1000*rand);


% TOGGLE
KF_ON = 1; 
CI_ON = 1; 
LE_ON = 1;
LKF_ON = 1;
CRLB_ON = 1;
TR_P_ON = 0;
FULL_ON = 0;


% PARAMS
m_vec = [1 2 3];

par = get_default_sim_params;
par.ID = 1;
par.sigma_r = 1000;
par.sigma_az = 1*d2r;
par.M = 10;
par.cntrl.globally_known_est = 1;


% --- SIMULATION ---
scen = load_scen(par); 
nagents = scen.nagents; ntargets = scen.ntargets;

itgt = 1; itgt = min([itgt ntargets]);
nm = length(m_vec);

if KF_ON; rgevkf = cell(nm,nagents); pgevkf = rgevkf; agevkf = rgevkf; end
if CI_ON; rgevci = cell(nm,nagents); pgevci = rgevci; agevci = rgevci; end
if LE_ON; rgevle = cell(nm,nagents); pgevle = rgevle; agevle = rgevle; end

if FULL_ON
    if KF_ON; rfkf = cell(1,nagents); pfkf = rfkf; afkf = rfkf; end
    if CI_ON; rfci = cell(1,nagents); pfci = rfci; afci = rfci; end
    if LE_ON; rfle = cell(1,nagents); pfle = rfle; afle = rfle; end
end

if LKF_ON; rlkf = cell(1,nagents); plkf = rlkf; alkf = rlkf; end

fprintf('\n\nSimulations starting at %s\n',get_datetime)

for im = 1:nm

    par.m = m_vec(im);
    fprintf('\n--- m = %d ---\n',par.m)

    % --- KF ---
    if KF_ON
        fprintf('\nRunning KF...\n')
        par.comm_mgmt_method = def_communication_management_technique('gevo');
        par.dimred_loss = def_dimred_loss_function('kf');
        par.fus_method = def_fusion_method('kf');
        rng(SEED)
        sim_out = simenv(par);
        
        for ia = 1:nagents
            [rgevkf{im,ia},pgevkf{im,ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
            [~,agevkf{im,ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
        end
    end
    
    % --- CI ---
    if CI_ON
        fprintf('\nRunning CI...\n')
        par.comm_mgmt_method = def_communication_management_technique('gevo');
        par.dimred_loss = def_dimred_loss_function('ci');
        par.fus_method = def_fusion_method('ci');
        rng(SEED) 
        sim_out = simenv(par);
        
        for ia = 1:nagents
            [rgevci{im,ia},pgevci{im,ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
            [~,agevci{im,ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
        end
    end
    
    % --- LE ---
    if LE_ON
        fprintf('\nRunning LE...\n')
        par.comm_mgmt_method = def_communication_management_technique('gevo');
        par.dimred_loss = def_dimred_loss_function('le');
        par.fus_method = def_fusion_method('le');
        rng(SEED)
        sim_out = simenv(par);
        
        for ia = 1:nagents
            [rgevle{im,ia},pgevle{im,ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
            [~,agevle{im,ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
        end
    end
end

if FULL_ON
    fprintf('\n--- m = p (full estimates) ---\n')
    if KF_ON 
        fprintf('\nRunning KF...\n')
        par.comm_mgmt_method = def_communication_management_technique('full');
        par.dimred_loss = def_dimred_loss_function('kf');
        par.fus_method = def_fusion_method('kf');
        rng(SEED)
        sim_out = simenv(par);

        for ia = 1:nagents
            [rfkf{ia},pfkf{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
            [~,afkf{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
        end
    end
    if CI_ON 
        fprintf('\nRunning CI...\n')
        par.comm_mgmt_method = def_communication_management_technique('full');
        par.dimred_loss = def_dimred_loss_function('ci');
        par.fus_method = def_fusion_method('ci');
        rng(SEED)
        sim_out = simenv(par);

        for ia = 1:nagents
            [rfci{ia},pfci{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
            [~,afci{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
        end
    end
    if LE_ON 
        fprintf('\nRunning LE...\n')
        par.comm_mgmt_method = def_communication_management_technique('full');
        par.dimred_loss = def_dimred_loss_function('le');
        par.fus_method = def_fusion_method('le');
        rng(SEED)
        sim_out = simenv(par);

        for ia = 1:nagents        
            [rfle{ia},pfle{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
            [~,afle{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
        end
    end
end

fprintf('\n--- MISC ---\n')

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

fprintf('\n')


% --- PLOT ---
nrow = 2;
ncol = nagents;

c = 1/par.sigma_r;

n = size(sim_out.rec{1}.xhat{1},1);
M = length(sim_out.rec{1}.P);
anees_ci = anees_confidence_interval(n,M);

N = scen.par.N; t = 1:N;
%N = 14;
idx = 1:1:N;
lw_vec = linspace(0.75,1.5,nm);
lw = 1.5;
lwb = 2.0;
lwp = 1.5;
lwf = 1.0;
ms = 10;
clr = get_colors_rgb;
clrkf = get_colormap('orange');
clrci = get_colormap('green'); %get_colormap('blue');
clrle = get_colormap('blue');
clrlkf = 0.4*[1 1 1];
cred = get_colormap('red');

figure(1);clf
for ia = 1:nagents

    subplot(nrow,ncol,ia); hold on

    if CRLB_ON; h = plot(t(2:end),c*pcrlb(2:end),'k-','DisplayName','CRLB'); h.LineWidth = lwb; end
    if LKF_ON; h = plot(t(idx),c*rlkf{ia}(idx),'k-','DisplayName','LKF'); h.Color = clrlkf; h.LineWidth = lwb; end
    
    if KF_ON
        for im = 1:nm
            h = plot(t(idx),c*rgevkf{im,ia}(idx),'-','DisplayName',sprintf('GEVO-KF $m = %d$',m_vec(im))); h.LineWidth = lw_vec(im); h.Color = clrkf(7,:);
            if TR_P_ON; h = plot(t(idx),c*pgevkf{im,ia}(idx),'--','DisplayName','GEVO-KF-P'); h.LineWidth = lwp; h.Color = clrkf(6,:); end
        end 
        if FULL_ON; h = plot(t(idx),c*rfkf{ia}(idx),'--','DisplayName','FULL-KF'); h.LineWidth = lwf; h.Color = clrkf(6,:); end
    end
    if CI_ON
        for im = 1:nm
            h = plot(t(idx),c*rgevci{im,ia}(idx),'-','DisplayName',sprintf('GEVO-CI $m = %d$',m_vec(im))); h.LineWidth = lw_vec(im); h.Color = clrci(7,:);
            if TR_P_ON; h = plot(t(idx),c*pgevci{im,ia}(idx),'--','DisplayName','GEVO-CI-P'); h.LineWidth = lwp; h.Color = clrci(6,:); end
        end 
        if FULL_ON; h = plot(t(idx),c*rfci{ia}(idx),'--','DisplayName','FULL-CI'); h.LineWidth = lwf; h.Color = clrci(6,:); end
    end
    if LE_ON
        for im = 1:nm
            h = plot(t(idx),c*rgevle{im,ia}(idx),'-','DisplayName',sprintf('GEVO-LE $m = %d$',m_vec(im))); h.LineWidth = lw_vec(im); h.Color = clrle(7,:);
            if TR_P_ON; h = plot(t(idx),c*pgevle{im,ia}(idx),'--','DisplayName','GEVO-LE-P'); h.LineWidth = lwp; h.Color = clrle(6,:); end
        end
        if FULL_ON; h = plot(t(idx),c*rfle{ia}(idx),'--','DisplayName','FULL-LE'); h.LineWidth = lwf; h.Color = clrle(6,:); end
    end

    ylim([0 c*rlkf{ia}(1)])
    xlabel('Time [s]','interpreter','latex'); 
    %set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    ylabel('RMSE$/\sigma_r$','interpreter','latex')
    title(sprintf('agent %d',ia))
    if ia == 1; legend show; end
    box on
end




% --- NEES ---
x = t(idx); x2 = [x fliplr(x)];

for ia = 1:nagents

    subplot(nrow,ncol,ia+ncol); 
    hold on

    %h = plot(t(idx),ones(1,length(idx)),'k-','DisplayName','ANEES = 1'); turn_legend_item_off(h);
    h = plot(t(idx),anees_ci.lvl99(1)*ones(1,length(idx)),'k--'); turn_legend_item_off(h);
    h = plot(t(idx),anees_ci.lvl99(2)*ones(1,length(idx)),'k--','DisplayName','99\% conf. int.'); 

    if LKF_ON; h = plot(t(idx),alkf{ia}(idx),'k-','DisplayName','LKF'); h.Color = clrlkf; h.LineWidth = lwb; end

    if KF_ON
        for im = 1:nm
            h = plot(t(idx),agevkf{im,ia}(idx),'-','DisplayName',sprintf('GEVO-KF $m = %d$',m_vec(im))); h.LineWidth = lw_vec(im); h.Color = clrkf(7,:);
        end
        if FULL_ON; h = plot(t(idx),afkf{ia}(idx),'--','DisplayName','FULL-KF'); h.LineWidth = lwf; h.Color = clrkf(6,:); end
    end

    if CI_ON
        for im = 1:nm
            h = plot(t(idx),agevci{im,ia}(idx),'-','DisplayName',sprintf('GEVO-CI $m = %d$',m_vec(im))); h.LineWidth = lw_vec(im); h.Color = clrci(7,:);
        end        
        if FULL_ON; h = plot(t(idx),afci{ia}(idx),'--','DisplayName','FULL-CI'); h.LineWidth = lwf; h.Color = clrci(6,:); end
    end

    if LE_ON
        for im = 1:nm
            h = plot(t(idx),agevle{im,ia}(idx),'-','DisplayName',sprintf('GEVO-LE $m = %d$',m_vec(im))); h.LineWidth = lw_vec(im); h.Color = clrle(7,:);
        end        
        if FULL_ON; h = plot(t(idx),afle{ia}(idx),'--','DisplayName','FULL-LE'); h.LineWidth = lwf; h.Color = clrle(6,:); end
    end

    %set(ax,'layer','top')
    ylim([0.25 2.25])
    xlabel('Time [s]','interpreter','latex'); ylabel('ANEES','interpreter','latex')
    if ia == 1; legend show; end
    box on

end




