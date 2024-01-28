function track_fusion_design
% --- track_fusion_design() -----------------------------------------------
% Example of evaluation of track fusion design. This example is used in the
% summary 
%
% "The Bright Side of Decentralized Target Tracking: A Thesis Summary 
% Dedicated to the Practitioner"
%
% Requirements: src/thesis_lib
%
% 2024-01-08 Robin Forsling

set_latex_interpreter;

SEED = round(1000*rand);


% TOGGLE
CI_ON = 1; 
LE_ON = 1;
NKF_ON = 1;
LKF_ON = 1;
CRLB_ON = 1;


% PARAMETERS
par = default_simulation_parameters;
par.ID = 2;
par.sigma_r = 1000;
par.sigma_az = 1.0*d2r;
par.M = 100;


% --- SIMULATIONS ---
scen = load_scen(par); 
nagents = scen.nagents; ntargets = scen.ntargets;

itgt = 1; itgt = min([itgt ntargets]);

s.rmse = cell(1,nagents);
s.rmt = cell(1,nagents);
s.anees = cell(1,nagents);
s.coin = cell(1,nagents);

CI = s;
LE = s;
NKF = s;
LKF = s;

fprintf('\n\nSimulations starting at %s\n',get_datetime)
fprintf('\n--- TRACK FUSION EVALUATION ---\n')

if CI_ON 
    fprintf('\nRunning CI...\n')
    par.comm_mgmt_method = def_communication_management_technique('full');
    par.fus_method = def_fusion_method('ci');
    rng(SEED)
    sim_out = simenv(par);

    for ia = 1:nagents
        [~,CI.rmt{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
        [~,~,CI.coin{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end
end

if LE_ON 
    fprintf('\nRunning LE...\n')
    par.comm_mgmt_method = def_communication_management_technique('full');
    par.fus_method = def_fusion_method('le');
    rng(SEED)
    sim_out = simenv(par);

    for ia = 1:nagents
        [~,LE.rmt{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
        [~,~,LE.coin{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end
end

if NKF_ON
    fprintf('\nRunning NKF...\n')    
    par.comm_mgmt_method = def_communication_management_technique('full');
    par.fus_method = def_fusion_method('kf');
    rng(SEED)
    sim_out = simenv(par);
    
    for ia = 1:nagents
        [~,NKF.rmt{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
        [~,~,NKF.coin{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
    end
end

if LKF_ON
    fprintf('\nRunning LKF...\n')    
    par.comm_mgmt_method = def_communication_management_technique('no-communication');
    par.dimred_loss = 0;
    par.fus_method = def_fusion_method('lkf');
    par.cntrl.globally_known_est = 0;  
    rng(SEED)
    sim_out = simenv(par);
    
    for ia = 1:nagents
        [~,LKF.rmt{ia}] = rmse(sim_out.scen.XT{itgt},sim_out.rec{ia});
        [~,~,LKF.coin{ia}] = nees(sim_out.scen.XT{itgt},sim_out.rec{ia});
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


%% EVALUATION

% --- PLOT ---
ia = 1;
N = scen.par.N; t = 1:N;
idx = 2:2:N;

nrow = 1;
ncol = 2;

n = size(sim_out.rec{1}.xhat{1},1);
M = length(sim_out.rec{1}.P);

lw = 1.5; lwb = 2.0; lwf = 2.0;
clr = get_thesis_colors;
clr.lkf = clr.darkgray;

if CRLB_ON; c = pcrlb(idx);
else; c = par.sigma_r;
end


figure(ia);clf


% --- RMT ---
subplot(nrow,ncol,1); hold on

if CI_ON; h = plot(t(idx),CI.rmt{ia}(idx)./c,'-','DisplayName','CI'); h.Color = clr.ci; h.LineWidth = lwb; end
if LE_ON; h = plot(t(idx),LE.rmt{ia}(idx)./c,'-','DisplayName','LE'); h.Color = clr.le; h.LineWidth = lwb; end
if NKF_ON; h = plot(t(idx),NKF.rmt{ia}(idx)./c,'-','DisplayName','NKF'); h.Color = clr.naive; h.LineWidth = lwb; end
if LKF_ON; h = plot(t(idx),LKF.rmt{ia}(idx)./c,'-','DisplayName','LKF'); h.Color = clr.lkf; h.LineWidth = lwb; end

ylim([0.0 4.0])
xlabel('$k$','interpreter','latex'); ylabel('RMT','interpreter','latex')
box on


% --- COIN ---
subplot(nrow,ncol,2); hold on

if CI_ON; h = plot(t(idx),CI.coin{ia}(idx),'-','DisplayName','CI'); h.Color = clr.ci; h.LineWidth = lwb; end
if LE_ON; h = plot(t(idx),LE.coin{ia}(idx),'-','DisplayName','LE'); h.Color = clr.le; h.LineWidth = lwb; end
if NKF_ON; h = plot(t(idx),NKF.coin{ia}(idx),'-','DisplayName','NKF'); h.Color = clr.naive; h.LineWidth = lwb; end
if LKF_ON; h = plot(t(idx),LKF.coin{ia}(idx),'-','DisplayName','LKF'); h.Color = clr.lkf; h.LineWidth = lwb; end

ylim([0.0 4.0])
xlabel('$k$','interpreter','latex'); ylabel('COIN','interpreter','latex')
legend show; 
box on

set_fontsize_all(14)



%% TIKZ

% --- RMT ---
fprintf('\\drawci\\plotcommand%s\n',get_tikz_plot_coordinates(t(idx),CI.rmt{ia}(idx)./c));
fprintf('\\drawle\\plotcommand%s\n',get_tikz_plot_coordinates(t(idx),LE.rmt{ia}(idx)./c));
fprintf('\\drawnkf\\plotcommand%s\n',get_tikz_plot_coordinates(t(idx),NKF.rmt{ia}(idx)./c));
fprintf('\\drawlkf\\plotcommand%s\n',get_tikz_plot_coordinates(t(idx),LKF.rmt{ia}(idx)./c));

fprintf('\n')

% --- COIN ---
fprintf('\\drawci\\plotcommand%s\n',get_tikz_plot_coordinates(t(idx),CI.coin{ia}(idx)));
fprintf('\\drawle\\plotcommand%s\n',get_tikz_plot_coordinates(t(idx),LE.coin{ia}(idx)));
fprintf('\\drawnkf\\plotcommand%s\n',get_tikz_plot_coordinates(t(idx),NKF.coin{ia}(idx)));
fprintf('\\drawlkf\\plotcommand%s\n',get_tikz_plot_coordinates(t(idx),LKF.coin{ia}(idx)));


