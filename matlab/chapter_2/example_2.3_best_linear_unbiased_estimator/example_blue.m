function example_blue
% --- example_blue() ------------------------------------------------------
% Example 2.3: Best Linear Unbiased Estimator
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- MODEL ---
H1 = eye(2); H2 = eye(2); H = [H1 ; H2];
R1 = [9 -2 ; -2 2];
R2 = [2 2 ; 2 9];
R12 = [1 1 ; -1 1];
R = [R1 R12 ; R12' R2];


% --- BLUE ---
Pblue = inv(H'/R*H);


% --- SAMPLED LINEAR ESTIMATORS ---
nsamp = 50;
Psamp = cell(nsamp,1);
for k = 1:nsamp
    S = get_random_covariance(4);
    K = inv(H'/S*H)*H'/S;
    Psamp{k} = K*R*K';
end


% --- PLOT ---
xc = [0;0];
clr = get_thesis_colors;
lw = 2;

figure(1);clf;hold on
for k = 1:nsamp; hsamp = plot_ellipse(xc,Psamp{k},'-'); hsamp.Color = clr.gray; end
h1 = plot_ellipse(xc,R1,'k-'); h1.Color = clr.darkyellow; h1.LineWidth = lw;
h2 = plot_ellipse(xc,R2,'k-'); h2.Color = clr.darkgreen; h2.LineWidth = lw;
hblue = plot_ellipse(xc,Pblue,'-'); hblue.Color = clr.blue; hblue.LineWidth = lw;

axis equal; box on
remove_ticks_and_ticklabels;
legend(gca,[h1 h2 hblue hsamp],'$R_1$','$R_2$','$P^\star$','$KRK^T,KH=I$')

set_fontsize_all(14)
