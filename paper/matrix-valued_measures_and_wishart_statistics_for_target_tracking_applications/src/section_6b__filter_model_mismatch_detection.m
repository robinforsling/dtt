function section_6b__filter_model_mismatch_detection
% --- section_6b__filter_model_mismatch_detection() -----------------------
% Application in Section 6-B: Filter model mismatch detection
%
% 2024-05-15 Robin Forsling


addpath('paper_lib/')
addpath(genpath('../../../src/thesis_lib/'))

set_latex_interpreter;


% --- SCENARIO ------------------------------------------------------------
par.M = 1e4;
par.nk = 60;
par.kvec = 1:par.nk;
par.kswitch = 10;

I2 = eye(2); O2 = zeros(2);


% ... PROCESS -------------------------------------------------------------
v0 = 10;
x0 = [0;0;v0;0];

qx = 1.5;
qy = 1/1.5;


% --- SENSOR --------------------------------------------------------------
H0 = [I2 O2];
R0 = 100*I2;


% --- ASSUMED PROCESS MODEL -----------------------------------------------
Ts = 1;
q = 10;

pm.Ts = Ts;
pm.q = q;
pm.F = [I2 Ts*I2 ; O2 I2];
pm.G = [Ts^2/2*I2 ; Ts*I2];
pm.Q = q^2*(pm.G*pm.G');


% --- ASSSUMED SENSOR MODEL -----------------------------------------------
sm.H = [I2 O2];
sm.R = R0; 


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
    y = generate_linear_measurement(xt,H0,R0);
    xhat = [y ; 0 ; 0];
    P = blkdiag(R0,v0^2*I2);

    for k = par.kvec
    

        % --- PROCESS ---
        q = pm.q; Ts = pm.Ts;
        C = diag([qx qy]);
        G = [Ts^2/2*C ; Ts*C];
        F = pm.F;

        w = q*G*randn(2,1);
        xt = F*xt + w;

        
        % --- MEASUREMENT ---
        y = generate_linear_measurement(xt,H0,R0);

    
        % --- TU ---
        xhat = pm.F*xhat;
        P = pm.F*P*pm.F' + pm.Q;


        % --- MU ---
        H = sm.H; R = sm.R;
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


% --- COMPUTE STATISTICS --------------------------------------------------
nis_stats = compute_nis_statistics(par,data);
Wishart_stats = compute_Wishart_statistics(par);
det_stats = compute_detections(par,nis_stats,Wishart_stats);


% --- RESULTS -------------------------------------------------------------
idx = 2:par.nk;
lw = 1;
clr = get_thesis_colors;
clrmax = clr.red; clrmean = clr.yellow; clrmin = clr.blue;

% --- NIS STATISTICS ---
figure(1);clf;hold on
plot_samples(par.kvec,nis_stats.lambdamean,clrmean);
plot_samples(par.kvec,nis_stats.lambdamin,clrmin);
plot_samples(par.kvec,nis_stats.lambdamax,clrmax);
h = plot(par.kvec(idx),Wishart_stats.lambdamin.cdf005(idx),'--'); h.Color = clrmin; h.LineWidth = lw;
h = plot(par.kvec(idx),Wishart_stats.lambdamax.cdf995(idx),'--'); h.Color = clrmax; h.LineWidth = lw;
h = plot(par.kvec(idx),Wishart_stats.lambdamean.cdf005(idx),'k--'); h.LineWidth = lw;
h = plot(par.kvec(idx),Wishart_stats.lambdamean.cdf995(idx),'k--'); h.LineWidth = lw;
hmin = plot(par.kvec,nis_stats.mean.lambdamin,'-'); hmin.Color = clrmin; hmin.LineWidth = lw;
hmax = plot(par.kvec,nis_stats.mean.lambdamax,'-'); hmax.Color = clrmax; hmax.LineWidth = lw;
hmean = plot(par.kvec,nis_stats.mean.lambdamean,'-'); hmean.Color = clrmean; hmean.LineWidth = lw;
legend(gca,[hmax hmean hmin],'$\lambda_{\max}$','$\bar{\lambda}$','$\lambda_{\min}$')
xlabel('$k$','interpreter','latex'); 

% --- DETECTIONS ---
figure(2);clf;hold on
hc = plot(par.kvec(idx),det_stats.chi2.p(idx),'-'); hc.Color = clrmean; hc.LineWidth = lw;
hw = plot(par.kvec(idx),det_stats.Wishart.p(idx),'-'); hw.Color = clrmax; hw.LineWidth = lw;
xlabel('$k$','interpreter','latex'); ylabel('$p_{out}$'); 

end




% --- LOCAL FUNCTIONS -----------------------------------------------------
function y = generate_linear_measurement(xt,H,R)
    L = chol(R,'lower');
    y = H*xt + L*randn(2,1);
end

function stats = compute_nis_statistics(par,data)
    m = 2;
    stats.M = par.M;
    stats.m = m;
    stats.nis = zeros(par.M,par.nk);
    stats.lambdamean = zeros(par.M,par.nk);
    stats.lambdamin = zeros(par.M,par.nk);
    stats.lambdamax = zeros(par.M,par.nk);
    for imc = 1:par.M
        Pi = zeros(m);
        for k = par.kvec
            ytilde = data{imc}.ytilde(:,k);
            L = chol(data{imc}.S(:,:,k),'lower'); 
            nu = L\ytilde;
            Pi = Pi + nu*nu';
            lambda = eig(Pi/k);
            stats.nis(imc,k) = sum(lambda);
            stats.lambdamean(imc,k) = stats.nis(imc,k)/m;
            stats.lambdamin(imc,k) = min(lambda);
            stats.lambdamax(imc,k) = max(lambda);
        end
    end
    stats.mean.nis = mean(stats.nis);
    stats.mean.lambdamean = mean(stats.lambdamean);
    stats.mean.lambdamin = mean(stats.lambdamin);
    stats.mean.lambdamax = mean(stats.lambdamax);
end

function stats = compute_Wishart_statistics(par)   
    m = 2; kvec = par.kvec; nk = length(kvec);
    stats.lambdamin.cdf005 = zeros(1,nk);
    stats.lambdamax.cdf995 = zeros(1,nk);
    stats.lambdamean.cdf005 = zeros(1,nk);
    stats.lambdamean.cdf995 = zeros(1,nk);
    for k = kvec
        n = k;
        if n >= m; stats.lambdamin.cdf005(k) = inv_cdf_lambdamin_Wishart(m,n,0.005)/n; end
        if n >= m; stats.lambdamax.cdf995(k) = inv_cdf_lambdamax_Wishart(m,n,0.995)/n; end
        stats.lambdamean.cdf005(k) = chi2inv(0.005,m*n)/(m*n);
        stats.lambdamean.cdf995(k) = chi2inv(0.995,m*n)/(m*n);
    end
end

function stats = compute_detections(par,nis_stats,Wishart_stats)
    kvec = par.kvec; nk = par.nk; M = par.M;
    cnt_Wishart = zeros(1,nk); cnt_chi2 = zeros(1,nk);
    for imc = 1:M
        for k = kvec
            if nis_stats.lambdamax(imc,k) >= Wishart_stats.lambdamax.cdf995(k) || nis_stats.lambdamin(imc,k) < Wishart_stats.lambdamin.cdf005(k)
                cnt_Wishart(k) = cnt_Wishart(k) + 1;
            end
            if nis_stats.lambdamean(imc,k) >= Wishart_stats.lambdamean.cdf995(k) || nis_stats.lambdamean(imc,k) < Wishart_stats.lambdamean.cdf005(k)
                cnt_chi2(k) = cnt_chi2(k) + 1;
            end
        end
    end
    stats.Wishart.cnt = cnt_Wishart;
    stats.Wishart.p = cnt_Wishart/M;
    stats.chi2.cnt = cnt_chi2;
    stats.chi2.p = cnt_chi2/M;
end

function plot_samples(kvec,lambda,clr)
    M = size(lambda,1);
    ns = min([M 10]);
    idx = 1:ns;
    lw = 0.5;
    opacity = 0.2;
    clr = [clr opacity];
    for i = 1:ns
        h = plot(kvec,lambda(idx(i),:),'-'); h.LineWidth = lw; h.Color = clr;
    end
end
