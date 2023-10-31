function varargout = example_pco_method
% --- example_pco_method() ------------------------------------------------
% Example 5.3: The Principal Component Optimization Method
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


Psi = [1 0];
sigmavec = 1:0.2:5;
ns = length(sigmavec);
tr_ratio = zeros(1,ns);

for i = 1:ns
    s = sigmavec(i);
    R1 = blkdiag(s^2,1);
    R2 = blkdiag(1,s^2);
    P = inv(inv(R1) + Psi'/(Psi*R2*Psi')*Psi);
    Pf = inv(inv(R1) + inv(R2));
    tr_ratio(i) = trace(P)/trace(Pf);
end

clr = get_thesis_colors;
lw = 1.5;

figure(1);clf
h = plot(sigmavec,tr_ratio,'-'); h.Color = clr.cyan; h.LineWidth = lw;
xlabel('$\sigma$','interpreter','latex')
ylabel('tr$(P)$/tr$(P_\mathrm{full})$','interpreter','latex')

set_fontsize_all(14)

if nargout > 0
    varargout{1} = get_tikz_plot_coordinates(sigmavec,tr_ratio);
end

