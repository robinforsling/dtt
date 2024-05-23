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

I2 = eye(2); O2 = zeros(2);


% --- PROCESS -------------------------------------------------------------
v0 = 10;
x0 = [0;0;v0;0];

Ts = 1;
q = 10;
alpha = 2;

Fc = [O2 I2 ; O2 O2];
Qc = q^2*[O2 O2 ; O2 diag([alpha 1/alpha])^2];
[F0,Q0] = process_continuous_to_discrete(Fc,Qc,Ts); 
Q0s = sqrtm(Q0);


% --- SENSOR --------------------------------------------------------------
H0 = [I2 O2];
R0 = 100*I2;


% --- ASSUMED PROCESS MODEL -----------------------------------------------
pm.Ts = Ts;
pm.q = q;
pm.F = [I2 Ts*I2 ; O2 I2];
pm.Q = q^2*[Ts^3*I2/3 Ts^2*I2/2 ; Ts^2*I2/2 Ts*I2];


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
        w = Q0s*randn(4,1);
        xt = F0*xt + w;

        
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
u_stats = compute_eigenvector_statistics(par,data);


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
xlabel('$k$','interpreter','latex'); ylabel('$p_{out}$','interpreter','latex'); 

% --- EIGENVECTORS ---
figure(3);clf;

subplot(1,2,1);hold on
plot_mean_and_std(par.kvec,u_stats.thetamin,clrmin)
xlabel('$k$','interpreter','latex'); ylabel('$\theta$ [degrees]','interpreter','latex')
title('$u_{\min}$')

subplot(1,2,2);hold on
plot_mean_and_std(par.kvec,u_stats.thetamax,clrmax)
xlabel('$k$','interpreter','latex'); ylabel('$\theta$ [degrees]','interpreter','latex')
title('$u_{\max}$')

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
        Pihat = zeros(m);
        for k = par.kvec
            ytilde = data{imc}.ytilde(:,k);
            B = chol(data{imc}.S(:,:,k),'lower'); 
            nu = B\ytilde;
            Pihat = Pihat + nu*nu';
            lambda = eig(Pihat/k);
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

function stats = compute_eigenvector_statistics(par,data)
    m = 2; kvec = par.kvec; nk = length(kvec);
    stats.mc.umin = zeros(2,nk);
    stats.mc.umax = zeros(2,nk);
    stats.run = cell(par.M);
    % Over MC:
    for k = kvec
        Pihat = zeros(2);
        for imc = 1:par.M
            ytilde = data{imc}.ytilde(:,k);
            B = chol(data{imc}.S(:,:,k),'lower');  
            nu = B\ytilde;
            Pihat = Pihat + nu*nu';
        end
        [U,Lambda] = eig(Pihat/par.M);
        [~,imin] = min(diag(Lambda));
        [~,imax] = max(diag(Lambda));
        umin = U(:,imin);
        umax = U(:,imax);
        stats.mc.umin(:,k) = B*umin/norm(B*umin);
        stats.mc.umax(:,k) = B*umax/norm(B*umax);
    end
    % Over single-run:
     theta.min = zeros(par.M,nk);
     theta.max = zeros(par.M,nk);
    for imc = 1:par.M
        stats.run{imc}.umin = zeros(2,nk);
        stats.run{imc}.umax = zeros(2,nk);
        Pihat = zeros(2);
        for k = kvec
            ytilde = data{imc}.ytilde(:,k);
            B = chol(data{imc}.S(:,:,k),'lower'); 
            nu = B\ytilde;
            Pihat = Pihat + nu*nu';
            [U,Lambda] = eig(Pihat/k);
            [~,imin] = min(diag(Lambda));
            [~,imax] = max(diag(Lambda));
            umin = U(:,imin); bmin = B*umin/norm(B*umin);
            umax = U(:,imax); bmax = B*umax/norm(B*umax);
            stats.run{imc}.umin(:,k) = bmin;
            stats.run{imc}.umax(:,k) = bmax;
            theta.min(imc,k) = angle_from_yaxis(bmin);
            theta.max(imc,k) = angle_from_xaxis(bmax);
        end
    end
    % Single-run theta over MC:
    stats.thetamin.mean = zeros(1,nk);
    stats.thetamin.std = zeros(1,nk);
    stats.thetamax.mean = zeros(1,nk);
    stats.thetamax.std = zeros(1,nk);
    for k = kvec
        stats.thetamin.mean(k) = mean(theta.min(:,k));
        stats.thetamin.std(k) = std(theta.min(:,k));
        stats.thetamax.mean(k) = mean(theta.max(:,k));
        stats.thetamax.std(k) = std(theta.max(:,k));
    end
end

function theta = angle_from_xaxis(v)
    v = sign(v(1))*v;
    theta = atan2(v(2),v(1));
end

function theta = angle_from_yaxis(v)
    v = sign(v(2))*v;
    theta = atan2(v(2),v(1))-pi/2;
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

function varargout = plot_mean_and_std(kvec,s,clr)
    lo = r2d*(s.mean-s.std); hi = r2d*(s.mean+s.std); av = r2d*s.mean;
    hconf = plot_confidence_interval(kvec,lo,hi,clr);
    hmean = plot(kvec,av,'-'); hmean.Color = clr; hmean.LineWidth = 2;
    if nargout > 0; varargout{1} = hmean; end
    if nargout > 1; varargout{2} = hconf; end
end

