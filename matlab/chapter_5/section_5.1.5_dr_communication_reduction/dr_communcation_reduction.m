function dr_communcation_reduction
% --- dr_communcation_reduction() -----------------------------------------
% Communication reduction when using dimension-reduced (DR) estimates.
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;

nfull = @(n2) n2.*(n2+3)/2;
ndr = @(m,n2) (2*m*n2-m^2+3*m)/2;
ndca = @(n2) 2*n2;
n2 = 2:20;
m = 1:4;


% --- PLOT ---
clr = get_thesis_colors;
lw = 1.5;

figure(1);clf;

subplot(1,2,1); hold on
for i = 1:length(m) 
    hdr = plot(n2,ndr(m(i),n2),'-','DisplayName',sprintf('$n_{\\mathrm{dr}}$, $m=%d$',m(i))); hdr.Color = clr.m{i}; hdr.LineWidth = lw;
end
hdca = plot(n2,ndca(n2),'k--','DisplayName','$n_{\mathrm{dca}}$'); hdca.LineWidth = lw;
hfull = plot(n2,nfull(n2),'-','DisplayName','$n_{\mathrm{full}}$'); hfull.Color = clr.darkgray; hfull.LineWidth = lw;
box on
xlabel('$n_2$','interpreter','latex')
ylabel('Number of parameters','interpreter','latex')
legend('show','interpreter','latex','location','northwest');

subplot(1,2,2); hold on
for i = 1:length(m)
    h = plot(n2,ndr(m(i),n2)./nfull(n2),'-','DisplayName',sprintf('$m=%d$',m(i))); h.Color = clr.m{i}; h.LineWidth = lw;
end
hdca = plot(n2,ndca(n2)./nfull(n2),'k--','DisplayName','$n_{\mathrm{dca}}$'); hdca.LineWidth = lw;
box on
xlabel('$n_2$','interpreter','latex')
ylabel('$n_{\mathrm{dr}}/n_{\mathrm{full}}$','interpreter','latex')
legend('show','interpreter','latex','location','northeast');

set_fontsize_all(14)





