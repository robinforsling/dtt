function example_common_information
% --- example_common_information() ----------------------------------------
% Example 3.8: Common Information
% 
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- SETTINGS ---
sigmae2 = 1;
sigmaw2_vec = [0.1 0.2 0.5 1 2 5 10];
nsw2 = length(sigmaw2_vec);

nk = 10;
k_vec = 1:nk;

rho_tu = zeros(nsw2,nk);
rho_mu = rho_tu; rho_fus = rho_tu;

for m = 1:nsw2
    
    R = diag(sigmae2*[1 1]);
    sigmaw2 = sigmaw2_vec(m);

    for k = 1:nk

        % TIME UPDATE:
        R = R + sigmaw2*ones(2);
        rho_tu(m,k) = compute_rho(R);

        % MEASUREMENT UPDATE:
        a = R(1,1)/(R(1,1)+sigmae2); b = R(2,2)/(R(2,2)+sigmae2);
        R = [a*R(1,1) a*b*R(1,2) ; a*b*R(1,2) b*R(2,2)];
        rho_mu(m,k) = compute_rho(R);

        % FUSION UPDATE:
        [P,R12] = bsc_fusion(R); 
        R(1,2) = R12; R(2,1) = R12;
        if mod(k,2) == 0; R(1,1) = P; 
        elseif mod(k,2) == 1; R(2,2) = P; 
        end
        rho_fus(m,k) = compute_rho(R);

    end
end


% --- PLOT ---
clr = get_thesis_colors;
lw = 2;

figure(1);clf;

subplot(1,3,1);hold on
for m = 1:nsw2
    q = sigmaw2_vec(m)/sigmae2;
    h = plot(k_vec,rho_tu(m,:),'-','DisplayName',sprintf('$q=%2.1f$',q)); h.Color = clr.spectral7{m}; h.LineWidth = lw;
end
ylim([0 1])
xlabel('$k$','interpreter','latex'); ylabel('$\rho_{k|k-1}$','interpreter','latex')
box on
legend('show','location','southeast'); 

subplot(1,3,2);hold on
for m = 1:nsw2
    h = plot(k_vec,rho_mu(m,:),'-'); h.Color = clr.spectral7{m}; h.LineWidth = lw;
end
ylim([0 1])
xlabel('$k$','interpreter','latex'); ylabel('$\rho_{k|k}$','interpreter','latex')
box on

subplot(1,3,3);hold on
for m = 1:nsw2
    h = plot(k_vec,rho_fus(m,:),'-'); h.Color = clr.spectral7{m}; h.LineWidth = lw;
end
ylim([0 1])
xlabel('$k$','interpreter','latex'); ylabel('$\rho_{k|k,f}$','interpreter','latex')
box on

set_fontsize_all(14)


end


% --- FUNCTIONS -----------------------------------------------------------
function rho = compute_rho(R)
    rho = R(1,2)/(sqrt(R(1,1))*sqrt(R(2,2)));
end

function [P,R12] = bsc_fusion(R)
    S = (R(1,1)+R(2,2)-2*R(1,2));
    K2 = (R(1,1)-R(1,2))/S;
    K1 = 1-K2;
    P = R(1,1) - K2*S*K2;
    R12 = K1*R(1,2)+K2*R(2,2);
end

