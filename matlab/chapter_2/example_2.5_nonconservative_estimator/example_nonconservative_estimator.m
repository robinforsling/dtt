function example_nonconservative_estimator
% --- example_nonconservative_estimator() ---------------------------------
% Example 2.5: Non-Conservative Estimator
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- MODEL ---
H1 = eye(2); H2 = eye(2); H = [H1 ; H2];
R1 = [9 -2 ; -2 2];
R2 = [2 2 ; 2 9];
R12 = [1 1 ; -1 1];
R = [R1 R12 ; R12' R2];
Rn = blkdiag(R1,R2);


% --- BLUE ---
Pblue = inv(H'/R*H);


% --- NAIVE ---
Pn = inv(H'/Rn*H);
Kn = inv(H'/Rn*H)*H'/Rn;
Pn0 = Kn*R*Kn';


% --- PLOT ---
xc = [0;0];
clr = get_thesis_colors;
lw = 2;

figure(1);clf;hold on
h1 = plot_ellipse(xc,R1,'k-'); h1.Color = clr.darkyellow; h1.LineWidth = lw;
h2 = plot_ellipse(xc,R2,'k-'); h2.Color = clr.darkgreen; h2.LineWidth = lw;
hblue = plot_ellipse(xc,Pblue,'-'); hblue.Color = clr.blue;
hn = plot_ellipse(xc,Pn,'-'); hn.Color = clr.naive; hn.LineWidth = lw;
hn0 = plot_ellipse(xc,Pn0,'--'); hn0.Color = clr.naive; hn0.LineWidth = lw;

axis equal; box on
remove_ticks_and_ticklabels;
legend(gca,[h1 h2 hblue hn hn0],'$R_1$','$R_2$','$P^\star$','$P$','cov$(\tilde{x})$')

set_fontsize_all(14)
