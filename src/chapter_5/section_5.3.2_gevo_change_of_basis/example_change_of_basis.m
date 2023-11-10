function example_change_of_basis
% --- example_change_of_basis() -------------------------------------------
% Section 5.3.2 The Generalized Eigenvalue Optimization Method - Change of 
% Basis Example: Figure 5.6
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;

R1 = [3 1 0 ; 1 2 0 ; 0 0 1];
R2 = [3 0 1 ; 0 2 0 ; 1 0 1];


% --- GEVO ---
m = 2;
Q = R1*R1; S = R1+R2;
[X,G] = eig(Q,S);
idx = get_max_idx_vec(G,m);
Phi = X(:,idx)';
Omega = gram_schmidt_process(Phi);

[U,~] = eig(Omega*R2*Omega');
Psi = U'*Omega;

RPhi = Phi*R2*Phi';
ROmega = Omega*R2*Omega';
RPsi = Psi*R2*Psi';

W = Psi';
IPhi = W'*Phi'/RPhi*Phi*W;
IOmega = W'*Omega'/ROmega*Omega*W;
IPsi = W'*Psi'/RPsi*Psi*W;


% --- PLOT ---
x0 = [0;0];
clr = get_thesis_colors;
lw = 2;

figure(1);clf;hold on
hRPhi = plot_ellipse(x0,RPhi,'-'); hRPhi.Color = clr.orange; hRPhi.LineWidth = lw;
hROmega = plot_ellipse(x0,ROmega,'-'); hROmega.Color = clr.darkyellow; hROmega.LineWidth = lw;
hRPsi = plot_ellipse(x0,RPsi,'-'); hRPsi.Color = clr.green; hRPsi.LineWidth = lw;
hIPhi = plot_ellipse(x0,IPhi,'k-'); hIPhi.LineWidth = 2.0;
hIOmega = plot_ellipse(x0,IOmega,'k-'); hIOmega.LineWidth = 1.5;
hIPsi = plot_ellipse(x0,IPsi,'k-'); hIPsi.LineWidth = 1.0;

%xlabel('$x_1$','interpreter','latex'); ylabel('$x_2$','interpreter','latex')
box on; axis equal
legend(gca,[hRPhi hROmega hRPsi],'$R_{\Phi}$','$R_{\Omega}$','$R_{\Psi}$')
remove_ticks_and_ticklabels;

set_fontsize_all(14)


