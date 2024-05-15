function section_5f__switching_target_dynamics
% --- section_5f__switching_target_dynamics() -----------------------------
% Switching target dynamics example of Section 5-F
%
% 2024-05-15 Robin Forsling


addpath('paper_lib/')
addpath(genpath('../../../src/thesis_lib/'))

set_latex_interpreter;


% --- SCENARIO ------------------------------------------------------------
par.M = 2e4; 
par.nk = 20;
par.kvec = 1:par.nk;
par.kswitch = 10;
Tk = 1;

T90 = [0 -1 ; 1 0]; % rotation matrix

I2 = eye(2); 
O2 = zeros(2);


% --- TARGET DYNAMICS -----------------------------------------------------
v0 = 10;
x0 = [0;0;v0;0];

q = 1;
F = [I2 Tk*I2 ; O2 I2];

% Before kswitch
G1 = [Tk^2/2*I2 ; Tk*I2];
Q1 = q^2*(G1*G1');

% From kswitch
qlong = sqrt(2);
qlat = 0.001;


% --- SENSOR --------------------------------------------------------------
H = [I2 O2];
R = 100*I2;


% --- MAIN LOOP -----------------------------------------------------------
data = cell(par.M,1);
for imc = 1:par.M

    mc.xt = zeros(4,par.nk);
    mc.xhat = mc.xt;
    mc.P = zeros(4,4,par.M);
    mc.ytilde = zeros(2,par.nk);
    mc.S = zeros(2,2,par.nk);

    xt = x0;


    % --- INITIALIZATION ---
    y = generate_linear_measurement(xt,H,R);
    xhat = [y ; 0 ; 0];
    P = blkdiag(R,v0^2*I2);
    for l = 1:4
        w = q*G1*randn(2,1);
        xt = F*xt + w;
        y = generate_linear_measurement(xt,H,R);
        xhat = F*xhat;
        P = F*P*F' + Q1;
        ytilde = y - H*xhat;
        S = H*P*H' + R;
        K = P*H'/S;
        xhat = xhat + K*ytilde;
        P = (eye(4)-K*H)*P;
    end


    for k = par.kvec
    

        % --- PROCESS SIMULATION ---
        if k < par.kswitch
            G = G1;
            Q = Q1;
        else
            ulong = xt(3:4)/norm(xt(3:4));
            ulat = T90*ulong;
            U = [ulong ulat];
            C = U*diag([qlat qlong]);
            G = [Tk^2/2*C ; Tk*C];
            Q = q^2*(G*G');
        end
        w = q*G*randn(2,1);
        xt = F*xt + w;

        
        % --- MEASUREMENT ---
        y = generate_linear_measurement(xt,H,R);

    
        % --- TU ---
        xhat = F*xhat;
        P = F*P*F' + Q1;


        % --- MU ---
        ytilde = y - H*xhat;
        S = H*P*H' + R;
        K = P*H'/S;
        xhat = xhat + K*ytilde;
        P = (eye(4)-K*H)*P;


        % --- SAVE DATA ---
        mc.xt(:,k) = xt;
        mc.xhat(:,k) = xhat;
        mc.P(:,:,k) = P;
        mc.ytilde(:,k) = ytilde;
        mc.S(:,:,k) = S;
        

    end

    data{imc} = mc;
end


% --- COMPUTE STATS -------------------------------------------------------
nees_stats = compute_nees_statistics(par,data);
nis_stats = compute_nis_statistics(par,data);
Wishart_stats = compute_Wishart_statistics(par,nees_stats);


% --- RESULTS -------------------------------------------------------------
lw = 1;
clr = get_thesis_colors;
clrmax = clr.red; clrmean = clr.yellow; clrmin = clr.blue;

figure(1);clf;

% NEES STATISTICS
subplot(1,2,1);hold on
plot_constant_line(par.kvec,Wishart_stats.nees.lambdamin.cdf005,clrmin,'--',lw); 
plot_constant_line(par.kvec,Wishart_stats.nees.lambdamax.cdf995,clrmax,'--',lw); 
plot_constant_line(par.kvec,Wishart_stats.nees.lambdamean.cdf005,clrmean,'--',lw); 
plot_constant_line(par.kvec,Wishart_stats.nees.lambdamean.cdf995,clrmean,'--',lw); 
hmin = plot(par.kvec,nees_stats.lambdamin,'-'); hmin.Color = clrmin; hmin.LineWidth = lw;
hmax = plot(par.kvec,nees_stats.lambdamax,'-'); hmax.Color = clrmax; hmax.LineWidth = lw;
hmean = plot(par.kvec,nees_stats.anees,'-'); hmean.Color = clrmean; hmean.LineWidth = lw;
legend(gca,[hmax hmean hmin],'$\lambda_{\max}$','$\bar{\lambda}$','$\lambda_{\min}$')
xlabel('$k$','interpreter','latex'); title('NEES Statistics')

% NIS STATISTICS
subplot(1,2,2);hold on
plot_constant_line(par.kvec,Wishart_stats.nis.lambdamin.cdf005,clrmin,'--',lw); 
plot_constant_line(par.kvec,Wishart_stats.nis.lambdamax.cdf995,clrmax,'--',lw); 
plot_constant_line(par.kvec,Wishart_stats.nis.lambdamean.cdf005,clrmean,'--',lw); 
plot_constant_line(par.kvec,Wishart_stats.nis.lambdamean.cdf995,clrmean,'--',lw); 
hmin = plot(par.kvec,nis_stats.lambdamin,'-'); hmin.Color = clrmin; hmin.LineWidth = lw;
hmax = plot(par.kvec,nis_stats.lambdamax,'-'); hmax.Color = clrmax; hmax.LineWidth = lw;
hmean = plot(par.kvec,nis_stats.anis,'-'); hmean.Color = clrmean; hmean.LineWidth = lw;
legend(gca,[hmax hmean hmin],'$\lambda_{\max}$','$\bar{\lambda}$','$\lambda_{\min}$')
xlabel('$k$','interpreter','latex'); title('NIS Statistics')

end




% --- LOCAL FUNCTIONS -----------------------------------------------------
function y = generate_linear_measurement(xt,H,R)
    L = chol(R,'lower');
    y = H*xt + L*randn(2,1);
end

function stats = compute_nees_statistics(par,data)
    idx = 1:4; 
    m = length(idx);
    stats.M = par.M;
    stats.m = m;
    stats.idx = idx;
    stats.Xi = zeros(m,m,par.nk);
    stats.nees = zeros(1,par.nk);
    stats.anees = zeros(1,par.nk);
    stats.lambdamin = zeros(1,par.nk);
    stats.lambdamax = zeros(1,par.nk);
    for k = par.kvec
        Xi = zeros(m);
        for imc = 1:par.M
            xtilde = data{imc}.xhat(idx,k) - data{imc}.xt(idx,k);
            L = chol(data{imc}.P(idx,idx,k),'lower'); 
            nu = L\xtilde;
            Xi = Xi + nu*nu';
        end
        Xi = Xi/par.M;
        lambda = eig(Xi);
        stats.Xi(:,:,k) = Xi;
        stats.nees(k) = trace(Xi);
        stats.anees(k) = stats.nees(k)/m;
        stats.lambdamin(k) = min(lambda);
        stats.lambdamax(k) = max(lambda);
    end
end

function stats = compute_nis_statistics(par,data)
    m = 2;
    stats.M = par.M;
    stats.m = m;
    stats.Pi = zeros(m,m,par.nk);
    stats.nis = zeros(1,par.nk);
    stats.anis = zeros(1,par.nk);
    stats.lambdamin = zeros(1,par.nk);
    stats.lambdamax = zeros(1,par.nk);
    for k = par.kvec
        Pi = zeros(m);
        for imc = 1:par.M
            ytilde = data{imc}.ytilde(:,k);
            L = chol(data{imc}.S(:,:,k),'lower'); 
            nu = L\ytilde;
            Pi = Pi + nu*nu';
        end
        Pi = Pi/par.M;
        lambda = eig(Pi);
        stats.Xi(:,:,k) = Pi;
        stats.nis(k) = trace(Pi);
        stats.anis(k) = stats.nis(k)/m;
        stats.lambdamin(k) = min(lambda);
        stats.lambdamax(k) = max(lambda);
    end
end

function stats = compute_Wishart_statistics(par,nees_stat)
    
    % NEES:
    m = nees_stat.m; n = par.M;
    stats.nees.lambdamin.cdf005 = inv_cdf_lambdamin_Wishart(m,n,0.005)/n;
    stats.nees.lambdamin.ev = ev_lambdamin_Wishart(m,n)/n;
    stats.nees.lambdamax.cdf995 = inv_cdf_lambdamax_Wishart(m,n,0.995)/n;
    stats.nees.lambdamax.ev = ev_lambdamax_Wishart(m,n)/n;
    stats.nees.lambdamean.cdf005 = chi2inv(0.005,m*n)/(m*n);
    stats.nees.lambdamean.cdf995 = chi2inv(0.995,m*n)/(m*n);

    % NIS:
    m = 2; n = par.M;
    stats.nis.lambdamin.cdf005 = inv_cdf_lambdamin_Wishart(m,n,0.005)/n;
    stats.nis.lambdamin.ev = ev_lambdamin_Wishart(m,n)/n;
    stats.nis.lambdamax.cdf995 = inv_cdf_lambdamax_Wishart(m,n,0.995)/n;
    stats.nis.lambdamax.ev = ev_lambdamax_Wishart(m,n)/n;
    stats.nis.lambdamean.cdf005 = chi2inv(0.005,m*n)/(m*n);
    stats.nis.lambdamean.cdf995 = chi2inv(0.995,m*n)/(m*n);

end

function varargout = plot_constant_line(kvec,val,clr,ls,lw)
    h = plot(kvec,val*ones(size(kvec)));
    h.Color = clr;
    h.LineStyle = ls;
    h.LineWidth = lw;
    if nargout > 0; varargout{1} = h; end
end

