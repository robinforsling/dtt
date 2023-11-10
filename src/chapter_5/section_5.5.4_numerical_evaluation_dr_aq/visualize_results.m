function prob_data = visualize_results(data)
% --- visualize_results() -------------------------------------------------
% Visualize results for the optimization strategy evaluation.
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- PARAMETERS ----------------------------------------------------------
M = data.sim.M;
N = data.scen.N;
n = data.scen.n;
nc = data.sim.nc;
c_vec = data.sim.c_vec;


% --- EMPIRICAL PROBABILITIES ---------------------------------------------
full = compute_probabilities(data.full,N);
red_fus = compute_probabilities(data.red_fus,N);
red_asso = compute_probabilities(data.red_asso,N);


% --- PLOT SETTINGS -------------------------------------------------------
clr = get_thesis_colors;
cfull = clr.ddarkgray;
cred_fus = clr.orange;
cred_asso = clr.cyan;
clr1 = clr.agent1;
clr2 = clr.agent2;
lw.mean = 2;
lw.std = 1;
ls.mean = '-';
ls.std = '--';


figure(1);clf;

% --- CORRECT ASSIGNMENT --------------------------------------------------
subplot(1,2,1);hold on
hfull = plot_mean_std(c_vec,full.ca,cfull,lw,ls);
hred_fus = plot_mean_std(c_vec,red_fus.ca,cred_fus,lw,ls);
hred_asso = plot_mean_std(c_vec,red_asso.ca,cred_asso,lw,ls);
xlabel('$c$','interpreter','latex'); ylabel('PCA');
title('Correct Assignment')
legend(gca,[hfull hred_fus hred_asso],'Full estimates','DR fusion optimal','DR association quality')
ylim([0 1.1])


% --- INCORRECT ASSIGNMENT ------------------------------------------------
subplot(1,2,2);hold on
hfull = plot_mean_std(c_vec,full.ia,cfull,lw,ls);
hred_fus = plot_mean_std(c_vec,red_fus.ia,cred_fus,lw,ls);
hred_asso = plot_mean_std(c_vec,red_asso.ia,cred_asso,lw,ls);
xlabel('$c$','interpreter','latex'); ylabel('PIA');
title('Incorrect Assignment')
ylim([0 1.1])

set_fontsize_all(14)


% --- TARGETS -------------------------------------------------------------
c = 0.1;
figure(2);clf;hold on
for i = 1:N

    % PLOT
    x = c*data.scen.X{i}(1:2);
    R1 = c^2*data.est.R1{i}(1:2,1:2);
    R2 = c^2*data.est.R2{i}(1:2,1:2);
    plot(x(1),x(2),'kx')
    h1 = plot_ellipse(x,R1); h1.Color = clr1;
    h2 = plot_ellipse(x,R2); h2.Color = clr2;

    legend(gca,[h1 h2],'$R_{1(i)}$','$R_{2(i)}$')

    box on; axis equal
    remove_ticks_and_ticklabels;

    set_fontsize_all(14)
end


% --- PACK RESULTS --------------------------------------------------------
prob_data.full = full;
prob_data.red_fus = red_fus;
prob_data.red_asso = red_asso;

end



% --- MISC FUNCTIONS ------------------------------------------------------
function p = compute_probabilities(s,N)
    p.ca.mean = s.ca.mean/N;
    p.ca.std_hi = (s.ca.mean+s.ca.std)/N;
    p.ca.std_lo = (s.ca.mean-s.ca.std)/N;
    p.ia.mean = s.ia.mean/N;
    p.ia.std_hi = (s.ia.mean+s.ia.std)/N;
    p.ia.std_lo = (s.ia.mean-s.ia.std)/N;
end

function hmean = plot_mean_std(c_vec,p,clr,lw,ls)
    hmean = plot(c_vec,p.mean,ls.mean); hmean.Color = clr; hmean.LineWidth = lw.mean;
    h = plot(c_vec,p.std_hi,ls.std); h.Color = clr; h.LineWidth = lw.std;
    h = plot(c_vec,p.std_lo,ls.std); h.Color = clr; h.LineWidth = lw.std;
end

function varargout = plot_ellipse(x0,S,varargin)
    L = chol(S,'lower');
    N = 1000;
    a = linspace(0,2*pi,N);
    x = L*[cos(a) ; sin(a)];
    h = plot(x0(1)+x(1,:),x0(2)+x(2,:),varargin{:});
    if nargout > 0; varargout{1} = h; end
end