function varargout = dca_communcation_reduction
% --- dca_communcation_reduction() ----------------------------------------
% Communication reduction when using the diagonal covariance approximation 
% (DCA).
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- SETUP ---
nfull = @(n) n.*(n+3)/2;
ndca = @(n) 2*n;
nx = 2:20;
lw = 1.5;


% --- PLOT ---
clr = get_thesis_colors;

figure(1);clf;

subplot(1,2,1); hold on
hdca = plot(nx,ndca(nx),'-'); hdca.Color = clr.orange; hdca.LineWidth = lw;
hfull = plot(nx,nfull(nx),'-'); hfull.Color = clr.black; hfull.LineWidth = lw;
box on
xlabel('$n_x$','interpreter','latex'); ylabel('Number of parameters','interpreter','latex')
legend(gca,[hfull hdca],'$n_{\mathrm{full}}$','$n_{\mathrm{dca}}$','location','northwest','interpreter','latex')

subplot(1,2,2); hold on
hdca = plot(nx,ndca(nx)./nfull(nx),'-'); hdca.Color = clr.orange; hdca.LineWidth = lw;
box on
xlabel('$n_x$','interpreter','latex'); ylabel('$n_{\mathrm{dca}}/n_{\mathrm{full}}$','interpreter','latex')

set_fontsize_all(14)


% --- OUTPUT ---
if nargout > 0
    n_str.full = get_tikz_plot_coordinates(nx,nfull(nx));
    n_str.dca = get_tikz_plot_coordinates(nx,ndca(nx));
    n_str.ratio = get_tikz_plot_coordinates(nx,ndca(nx)./nfull(nx));
    varargout{1} = n_str;
end



