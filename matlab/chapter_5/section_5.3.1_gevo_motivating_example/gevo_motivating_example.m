function varargout = gevo_motivating_example
% --- gevo_motivating_example() -------------------------------------------
% Section 5.3.1 Motivating Example - GEVO: Figure 5.5
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;

R1 = [3.2 1.2 ; 1.2 1.8]; I1 = inv(R1);
R2 = [4 0 ; 0 1];
H = @(a) [cos(a) sin(a)];
avec = 0:5:180;
na = length(avec);

trP = zeros(1,na);
trPfull = trace(inv(I1+inv(R2)));


% --- COMPUTE TRACE P -----------------------------------------------------
for i = 1:na
    Psi = H(d2r*avec(i));
    RPsi = Psi*R2*Psi';
    P = inv(I1 + Psi'/RPsi*Psi);
    trP(i) = trace(P);
end


% --- PLOT ----------------------------------------------------------------
clr = get_thesis_colors;
lw = 1.5; 

figure(1);clf;hold on
h = plot(avec,trPfull*ones(1,na),'k-'); h.LineWidth = lw;
h = plot(avec,trP,'-'); h.Color = clr.cyan; h.LineWidth = lw;
box on
xlabel('$\alpha$ [deg]','interpreter','latex')
ylabel('tr$(P)$','interpreter','latex')

set_fontsize_all(14)

if nargout > 0; varargout{1} = get_tikz_plot_coordinates(avec,trP); end
if nargout > 1; s.max = max(trP); s.min = min(trP); varargout{2} = s; end