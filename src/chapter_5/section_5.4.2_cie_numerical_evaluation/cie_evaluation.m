function cie_evaluation()
% --- cie_evaluation() ----------------------------------------------------
% Section 5.4.2 Numerical Evaluation - CIE: Figure 5.15
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- PARAMETERS ----------------------------------------------------------
m_vec = [1 2];

par = default_simulation_parameters;
par.ID = 1;
par.sigma_r = 1000;
par.sigma_az = 1.0*d2r;
par.M = 100;
par.SEED = round(10000*rand);

itgt = 1;


% --- RUN COMPARISON ------------------------------------------------------
[data_par,data,crlb] = run_cie_comparison(par,m_vec);

scen = load_scen(par);
fprintf('\n--- AGENT PARAMETERS ---\n')
fprintf('r\t sr\t saz\t r*saz/sr\n')
for i = 1:length(scen.sensor_pos)
    R = scen.agents{i}.sensor_model.R;
    sr = sqrt(R(1,1));
    saz = sqrt(R(2,2));
    xrel = scen.targets{itgt}.X(1:2,1)-scen.sensor_pos{i};
    ri = norm(xrel);
    fprintf('%4.0f\t %d\t %2.4f\t %2.5f \n',ri,sr,saz,ri*saz/sr)
end
fprintf('\n')


% --- PLOT ----------------------------------------------------------------
ia = 1;
nm = length(m_vec);

ls_vec = {'-','--',':','-.','-','--',':','-.'};
lw = 1.5; lwb = 2.0; 
lwdca = 1;
clr = get_thesis_colors;
clrdkf = clr.darkyellow;
clrde = clr.orange;

nx = data_par.nx; M = data_par.M;
anees_ci = anees_confidence_interval(nx,M);

N = data_par.N; t = 0:N-1; t = t+1;
idx = 2:2:N;
if ri*saz > sr; c = 1/(ri*saz); else; c = 1/sr; end
%c = 1/par.sigma_r;
xlimits = [1 19];

res = data{ia};


% --- FIGURE ---
figure(1);clf


% RMTR-LTG
subplot(2,2,1); hold on

for im = 1:nm
    h = plot(t(idx),res.kf{im}.loc.rmt.pos(idx)./res.kf{im}.glob.rmt.pos(idx),'-','DisplayName',sprintf('dKF $m = %d$',m_vec(im))); set_h(h,clrdkf,ls_vec{im},lw);
    h = plot(t(idx),res.ci{im}.loc.rmt.pos(idx)./res.ci{im}.glob.rmt.pos(idx),'-','DisplayName',sprintf('CI $m = %d$',m_vec(im))); set_h(h,clr.ci,ls_vec{im},lw);
    h = plot(t(idx),res.le{im}.loc.rmt.pos(idx)./res.le{im}.glob.rmt.pos(idx),'-','DisplayName',sprintf('LE $m = %d$',m_vec(im))); set_h(h,clr.le,ls_vec{im},lw);
end 

ylim([0.95 1.1]); xlim(xlimits)
xlabel('$k$','interpreter','latex'); ylabel('RMTR-LTG','interpreter','latex')
legend('show'); box on


% ANEES
subplot(2,2,2); hold on

h = plot(t(idx),anees_ci.lvl99(1)*ones(1,length(idx)),'k--'); turn_legend_item_off(h);
hconfint = plot(t(idx),anees_ci.lvl99(2)*ones(1,length(idx)),'k--'); 
h = plot(t(idx),res.lkf.anees(idx),'k-'); h.Color = clr.lkf; h.LineWidth = lwb;

for im = 1:nm
    h = plot(t(idx),res.kf{im}.loc.anees(idx),'-'); set_h(h,clrdkf,ls_vec{im},lw);
    h = plot(t(idx),res.kf{im}.loc.anees(idx),'-'); set_h(h,clr.ci,ls_vec{im},lw);
    h = plot(t(idx),res.kf{im}.loc.anees(idx),'-'); set_h(h,clr.le,ls_vec{im},lw);
end 
h = plot(t(idx),res.dca.anees(idx),'-','DisplayName','DCA-EIG'); set_h(h,clrde,'--',lwdca);

ylim([0 2]); xlim(xlimits)
xlabel('$k$','interpreter','latex'); ylabel('ANEES','interpreter','latex')
legend(gca,[hconfint],'99\% conf. int.'); 
box on


% RMT
subplot(2,2,3); hold on

h = plot(t(idx),c*crlb.glob.rt.pos(idx),'k-','DisplayName','CRLB'); h.LineWidth = lwb;
h = plot(t(idx),c*res.lkf.rmt.pos(idx),'k-','DisplayName','LKF'); h.Color = clr.lkf; h.LineWidth = lwb; 
for im = 1:nm
    h = plot(t(idx),c*res.kf{im}.loc.rmt.pos(idx),'-','DisplayName',sprintf('dKF $m = %d$',m_vec(im))); set_h(h,clrdkf,ls_vec{im},lw);
    h = plot(t(idx),c*res.ci{im}.loc.rmt.pos(idx),'-','DisplayName',sprintf('CI $m = %d$',m_vec(im))); set_h(h,clr.ci,ls_vec{im},lw);
    h = plot(t(idx),c*res.le{im}.loc.rmt.pos(idx),'-','DisplayName',sprintf('LE $m = %d$',m_vec(im))); set_h(h,clr.le,ls_vec{im},lw);
end 
h = plot(t(idx),c*res.dca.rmt.pos(idx),'--','DisplayName','DCA-EIG'); set_h(h,clrde,'--',lwdca);

ylim([0.0 1.0]); xlim(xlimits)
xlabel('$k$','interpreter','latex'); ylabel('RMT/$\sigma$','interpreter','latex')
box on


% RMSE
subplot(2,2,4); hold on

h = plot(t(idx),c*crlb.glob.rt.pos(idx),'k-','DisplayName','CRLB'); h.LineWidth = lwb;
h = plot(t(idx),c*res.lkf.rmse.pos(idx),'k-','DisplayName','LKF'); h.Color = clr.lkf; h.LineWidth = lwb; 
for im = 1:nm
    h = plot(t(idx),c*res.kf{im}.loc.rmse.pos(idx),'-','DisplayName',sprintf('dKF $m = %d$',m_vec(im))); set_h(h,clrdkf,ls_vec{im},lw);
    h = plot(t(idx),c*res.ci{im}.loc.rmse.pos(idx),'-','DisplayName',sprintf('CI $m = %d$',m_vec(im))); set_h(h,clr.ci,ls_vec{im},lw);
    h = plot(t(idx),c*res.le{im}.loc.rmse.pos(idx),'-','DisplayName',sprintf('LE $m = %d$',m_vec(im))); set_h(h,clr.le,ls_vec{im},lw);
end 
h = plot(t(idx),c*res.dca.rmse.pos(idx),'-','DisplayName','DCA-EIG'); set_h(h,clrde,'--',lwdca);

ylim([0.0 1.0]); xlim(xlimits)
xlabel('$k$','interpreter','latex'); ylabel('RMSE/$\sigma$','interpreter','latex')
box on

set_fontsize_all(14)


end




% --- FUNCTIONS -----------------------------------------------------------
function set_h(h,clr,ls,lw)
    if ~isempty(clr); h.Color = clr; end
    if ~isempty(ls); h.LineStyle = ls; end
    if ~isempty(lw); h.LineWidth = lw; end
end

