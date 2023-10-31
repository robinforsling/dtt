function example_hyperrectangle_enclosing
% --- example_hyperrectangle_enclosing() ----------------------------------
% Section 3.3.2 Methods for Preserving Conservativeness - DCA: Figure 3.8
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- SETTINGS ---
d = [4 1]; 
a = sqrt(d(1)); b = sqrt(d(2));

w_vec = 0.2:0.1:0.5;
nw = length(w_vec);
D = cell(nw,1);

nsamp = 100;

w1 = zeros(nw,1); w2 = w1;
for i = 1:nw
    w = [w_vec(i) 1-w_vec(i)];
    w1(i) = w(1); w2(i) = w(2);
    D{i} = diag(d./w);
end


% --- PLOT ---
clr = get_thesis_colors;
lw = 2;
x0 = [0;0];

figure(1);clf;hold on
for i = 1:nsamp
    rho = sample_rho([-2 2]);
    S = diag(d); S(1,2) = rho; S(2,1) = rho;
    h = plot_ellipse(x0,S); h.Color = clr.gray;
    turn_legend_item_off(h);
end

for i = 1:nw
    h = plot_ellipse(x0,D{i},'DisplayName',sprintf('$$\\omega_1=%1.1f,\\omega_2=%1.1f$$',w1(i),w2(i)));
    h.Color = clr.m{i}; h.LineWidth = lw;
end
rectangle('Position',[-a -b 2*a 2*b])

axis equal; box on
remove_ticks_and_ticklabels;
legend show

set_fontsize_all(14)

end 


function rho = sample_rho(interval)
    L = interval(2)-interval(1);
    rho = L*rand + interval(1);
end

