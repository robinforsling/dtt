function example_conservative_and_nonconservative_estimator
% --- example_conservative_and_nonconservative_estimator() ----------------
% Example 4.3: Conservative and Non-Conservative Estimators
% 
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% SCENARIO
H = [eye(2);eye(2)];
R1 = [13 3 ; 3 3];
R2 = [13 -6 ; -6 6];
x0 = [0;0];


% COMPUTATIONS
R0 = [R1 0.5*magic(2)' ; 0.5*magic(2) R2];
C = cell(5,1); A = C;
C{1} = 2*eye(2);
C{2} = [3 1 ; 2 2];
C{3} = [4 -1 ; -2 1];
C{4} = [1 3 ; 2 1];
C{5} = [0.5 -1 ; -1 0.5];
for i = 1:length(C)
    R12 = C{i};
    R = [R1 R12 ; R12' R2];
    if min(eig(R)) <= 0; error('not good...')
    else; A{i} = R;
    end
end

I1 = inv(R1); I2 = inv(R2);
PCI = inv(0.5*I1+0.5*I2);
KCI = PCI*(0.5*[I1 I2]);

R = R0;
P = inv(H'/R*H);
K = P*H'/R;


% --- PLOT ---
clr = get_thesis_colors;
lw = 2;

figure(1);clf;

% CONSERVATIVE ESTIMATOR
subplot(1,2,1);hold on
for i = 1:length(A)
    R = A{i};
    P = KCI*R*KCI'; 
    h = plot_ellipse(x0,P,'--'); h.Color = clr.darkgray;
end
hci = plot_ellipse(x0,PCI,'k-'); hci.Color = clr.green; hci.LineWidth = lw;

box on; axis equal
title('Conservative estimator')
legend(gca,[h hci],'$KSK^T,S\in\mathcal{A}$','$P$','location','northeast')
remove_ticks_and_ticklabels;

% NON-CONSERVATIVE ESTIMATOR
subplot(1,2,2);hold on
for i = 1:length(A)
    R = A{i};
    P = K*R*K'; 
    h = plot_ellipse(x0,P,'--'); h.Color = clr.darkgray;
end
hnc = plot_ellipse(x0,P,'k-'); hnc.Color = clr.red; hnc.LineWidth = lw;

box on; axis equal
title('Non-conservative estimator')
legend(gca,[h hnc],'$K^nS(K^n)^T,S\in\mathcal{A}$','$P^n$','location','northeast')
remove_ticks_and_ticklabels;

set_fontsize_all(14)
