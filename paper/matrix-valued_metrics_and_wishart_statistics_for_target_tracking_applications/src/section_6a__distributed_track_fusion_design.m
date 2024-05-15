function section_6a__distributed_track_fusion_design
% --- section_6a__distributed_track_fusion_design() -----------------------
% Application in Section 6-A: Distributed track fusion design
%
% 2024-05-15 Robin Forsling


addpath('paper_lib/')
addpath(genpath('../../../src/thesis_lib/'))

set_latex_interpreter;


% --- SIMULATION PARAMETERS -----------------------------------------------
par.M = 1e4;
par.nk = 60;
par.kvec = 1:par.nk;
par.nagents = 2;

iagent_eval = 1;


% --- PROCESS -------------------------------------------------------------
xt0 = [2.0e3;10e3;-100;100];
xt = get_central_motion(1,xt0,par.nk+1);
xt = xt(:,2:end);

T = 1;
q = 5;
pm.T = T;
pm.F = [eye(2) T*eye(2) ; zeros(2) eye(2)];
pm.Q = q^2*[T^3/3*eye(2) T^2/2*eye(2) ; T^2/2*eye(2) T*eye(2)]; 
pm.v0 = 100;


% --- SENSOR --------------------------------------------------------------
sigmar = 100;
sigmaaz = 2*d2r;

sm.h = @(x) [sqrt(x(1)^2+x(2)^2) ; atan2(x(2),x(1))];
sm.J = @(x) [x(1)/sqrt(x(1)^2+x(2)^2) x(2)/sqrt(x(1)^2+x(2)^2) zeros(1,2) ; ...
            -x(2)/(x(1)^2+x(2)^2) x(1)/(x(1)^2+x(2)^2) zeros(1,2)];
 
sm.L = diag([sigmar sigmaaz]);
sm.R = sm.L^2;

xs = cell(par.nagents,1);
xs{1} = [0;0];
xs{2} = [5e3;0];


% --- MC LOOP -------------------------------------------------------------
data = get_record_struct(par);
for imc = 1:par.M


    dl_schedule = get_dl_schedule(par);


    % --- INITIALIZE ---
    agent = cell(par.nagents,1);
    for i = 1:par.nagents
        y = generate_measurement(sm,xs{i},xt0);
        agent{i} = initialization(pm,sm,xs{i},y);
    end
    

    for k = 1:par.nk    
        

        % --- FILTER ---
        for i = 1:par.nagents
            a = agent{i};
            y = generate_measurement(sm,xs{i},xt(:,k));
            [a.lkf,a.nkf,a.ci,a.le] = time_update(pm,a.lkf,a.nkf,a.ci,a.le);
            [a.lkf,a.nkf,a.ci,a.le] = measurement_update(sm,xs{i},y,a.lkf,a.nkf,a.ci,a.le);
            agent{i} = a;
        end


        % --- INTERACTION ---
        [idx_tx,idx_rx,data_tx] = simulate_interactions(agent,dl_schedule,k);
  

        % --- FUSION ---
        for itx = 1:length(idx_tx)
            est_tx = data_tx{itx};
            for irx = idx_rx{itx}
                est_rx = agent{irx};
                est_fus = data_fusion(est_rx,est_tx);
                agent{irx} = est_fus;
            end
        end


        % --- SAVE DATA ---
        for i = 1:par.nagents
            data.lkf{imc,i}.xhat(:,k) = agent{i}.lkf.xhat;
            data.lkf{imc,i}.P(:,:,k) = agent{i}.lkf.P;
            data.nkf{imc,i}.xhat(:,k) = agent{i}.nkf.xhat;
            data.nkf{imc,i}.P(:,:,k) = agent{i}.nkf.P;
            data.ci{imc,i}.xhat(:,k) = agent{i}.ci.xhat;
            data.ci{imc,i}.P(:,:,k) = agent{i}.ci.P;
            data.le{imc,i}.xhat(:,k) = agent{i}.le.xhat;
            data.le{imc,i}.P(:,:,k) = agent{i}.le.P;
        end
    end
end


% --- COMPUTE STATISTICS --------------------------------------------------
crlb = compute_crlb(par,sm,pm,xs,xt0,xt);
perf_stats = compute_performance_statistics(par,data,xt,iagent_eval);
nees_stats = compute_nees_statistics(par,data,xt,iagent_eval);
Wishart_stats = compute_Wishart_statistics(nees_stats.m,par.M,0.9995);



% --- RESULTS -------------------------------------------------------------
switch iagent_eval
    case 1; kidx = 2:2:par.nk;
    case 2; kidx = 1:2:par.nk;
end
clr = get_thesis_colors;


% --- PERFORMANCE ---
figure(1);clf;
plot_performance_stats(kidx,par,perf_stats,crlb,clr)


% --- CONSERVATIVENESS ---
figure(2);clf

subplot(1,3,1);hold on
hr = plot_nees_region(kidx,par,nees_stats.lkf,clr.lkf);
hlkf = plot_nees_statistics(kidx,par,nees_stats.lkf,clr.lkf);
hr = plot_nees_region(kidx,par,nees_stats.nkf,clr.naive);
hnkf = plot_nees_statistics(kidx,par,nees_stats.nkf,clr.naive);
hw = plot_Wishart_statistics(kidx,par,Wishart_stats,clr.black);
add_plot_description('NKF')

subplot(1,3,2);hold on
hr = plot_nees_region(kidx,par,nees_stats.lkf,clr.lkf);
hlkf = plot_nees_statistics(kidx,par,nees_stats.lkf,clr.lkf);
hr = plot_nees_region(kidx,par,nees_stats.ci,clr.ci);
hci = plot_nees_statistics(kidx,par,nees_stats.ci,clr.ci);
hw = plot_Wishart_statistics(kidx,par,Wishart_stats,clr.black);
add_plot_description('CI')

subplot(1,3,3);hold on
hr = plot_nees_region(kidx,par,nees_stats.lkf,clr.lkf);
hlkf = plot_nees_statistics(kidx,par,nees_stats.lkf,clr.lkf);
hr = plot_nees_region(kidx,par,nees_stats.le,clr.le);
hle = plot_nees_statistics(kidx,par,nees_stats.le,clr.le);
hw = plot_Wishart_statistics(kidx,par,Wishart_stats,clr.black);
add_plot_description('LE')

end




% --- TRACKING ------------------------------------------------------------
function y = generate_measurement(sm,xs,xt)
    y = sm.h(xt(1:2)-xs(1:2)) + sm.L*randn(2,1);
end

function agent = initialization(pm,sm,xs,y)
    r = y(1); az = y(2);
    x0 = [r*cos(az)+xs(1) ; r*sin(az)+xs(2) ; zeros(2,1)];
    J = [cos(az) -r*sin(az) ; sin(az) r*cos(az)];
    P0 = 1.5*blkdiag(J*sm.R*J',pm.v0^2*eye(2));
    est = get_estimate_struct(x0,P0);
    agent.lkf = est;
    agent.nkf = est;
    agent.ci = est;
    agent.le = est;
end

function varargout = time_update(pm,varargin)
    F = pm.F; Q = pm.Q;
    for i = 1:nargout
        est = varargin{i}; 
        est.xhat = F*est.xhat;
        est.P = F*est.P*F' + Q;
        varargout{i} = est;
    end 
end

function varargout = measurement_update(sm,xs,y,varargin)
    R = sm.R; 
    for i = 1:nargout
        xhat = varargin{i}.xhat; P = varargin{i}.P;
        xrel = xhat(1:2)-xs;
        H = sm.J(xrel);
        ytilde = y-sm.h(xrel);
        S = H*P*H' + R;
        K = P*H'/S;
        xhat = xhat + K*ytilde;
        P = (eye(4)-K*H)*P;
        varargout{i}.xhat = xhat;
        varargout{i}.P = P;
    end
end

function est_fus = data_fusion(est_rx,est_tx)

    est_fus = est_rx;
    

    % --- LKF (no fusion) ---


    % --- NKF ---
    [y1,R1] = extract_estimate(est_rx.nkf); [y2,R2] = extract_estimate(est_tx.nkf);
    [xhat,P] = naive_fusion(y1,R1,y2,R2);
    est_fus.nkf.xhat = xhat; est_fus.nkf.P = P; 


    % --- CI ---
    [y1,R1] = extract_estimate(est_rx.ci); [y2,R2] = extract_estimate(est_tx.ci);
    [xhat,P] = ci_fusion(y1,R1,y2,R2);
    est_fus.ci.xhat = xhat; est_fus.ci.P = P;


    % --- LE ---
    [y1,R1] = extract_estimate(est_rx.le); [y2,R2] = extract_estimate(est_tx.le);
    [xhat,P] = le_fusion(y1,R1,y2,R2);
    est_fus.le.xhat = xhat; est_fus.le.P = P;

end

function [xhat,P] = naive_fusion(y1,R1,y2,R2)
    I1 = inv(R1); I2 = inv(R2);
    P = inv(I1+I2);
    xhat = P*(I1*y1+I2*y2);
end

function [xhat,P] = ci_fusion(y1,R1,y2,R2)
    I1 = inv(R1); I2 = inv(R2);
    J = @(w) trace(inv(w*I1+(1-w)*I2));
    w = fminbnd(J,0,1);
    P = inv(w*I1+(1-w)*I2);
    xhat = P*(w*I1*y1+(1-w)*I2*y2);
end

function [xhat,P] = le_fusion(y1,R1,y2,R2)
    nx = length(y1);
    [U1,D1] = eig(R1); T1 = sqrtm(inv(D1))*U1';
    [U2,D2] = eig(T1*R2*T1'); T = U2'*T1;
    z1 = T*y1; z2 = T*y2;
    z = z1; D = eye(nx);
    for i = 1:nx
        if D2(i,i) <= 1; z(i) = z2(i); D(i,i) = D2(i,i); end
    end
    Ti = inv(T);
    xhat = Ti*z;
    P = Ti*D*Ti';
end

function [idx_tx,idx_rx,data_tx] = simulate_interactions(agent,dl_schedule,k)
    dl = dl_schedule(:,k);
    nagents = length(agent);
    idx_tx = get_nonzero_indices(dl);
    ntx = length(idx_tx);
    idx_rx = cell(1,ntx); 
    data_tx = cell(1,ntx);
    for i = 1:ntx
        idx_rx{i} = get_index_vector_diff(nagents,idx_tx{i});
        data_tx{i} = agent{idx_tx{i}};
    end
end

function [xhat,P] = extract_estimate(est)
    xhat = est.xhat; P = est.P;
end




% --- EVALUATION ----------------------------------------------------------
function stats = compute_performance_statistics(par,data,xt,iagent)

    s.rmt = zeros(1,par.nk);
    s.rmse = zeros(1,par.nk);
    stats.lkf = s;
    stats.nkf = s;
    stats.ci = s;
    stats.le = s;

    % --- LKF ---
    for k = par.kvec
        t = 0; e = 0;
        for imc = 1:par.M
            xtilde = data.lkf{imc,iagent}.xhat(:,k) - xt(:,k);
            P = data.lkf{imc,iagent}.P(:,:,k);
            e = e + xtilde'*xtilde;
            t = t + trace(P);
        end
        stats.lkf.rmse(k) = sqrt(e/par.M);
        stats.lkf.rmt(k) = sqrt(t/par.M);
    end

    % --- NKF ---
    for k = par.kvec
        t = 0; e = 0;
        for imc = 1:par.M
            xtilde = data.nkf{imc,iagent}.xhat(:,k) - xt(:,k);
            P = data.nkf{imc,iagent}.P(:,:,k);
            e = e + xtilde'*xtilde;
            t = t + trace(P);
        end
        stats.nkf.rmse(k) = sqrt(e/par.M);
        stats.nkf.rmt(k) = sqrt(t/par.M);
    end  

    % --- CI ---
    for k = par.kvec
        t = 0; e = 0;
        for imc = 1:par.M
            xtilde = data.ci{imc,iagent}.xhat(:,k) - xt(:,k);
            P = data.ci{imc,iagent}.P(:,:,k);
            e = e + xtilde'*xtilde;
            t = t + trace(P);
        end
        stats.ci.rmse(k) = sqrt(e/par.M);
        stats.ci.rmt(k) = sqrt(t/par.M);
    end 

    % --- LE ---
    for k = par.kvec
        t = 0; e = 0;
        for imc = 1:par.M
            xtilde = data.le{imc,iagent}.xhat(:,k) - xt(:,k);
            P = data.le{imc,iagent}.P(:,:,k);
            e = e + xtilde'*xtilde;
            t = t + trace(P);
        end
        stats.le.rmse(k) = sqrt(e/par.M);
        stats.le.rmt(k) = sqrt(t/par.M);
    end 
end

function crlb = compute_crlb(par,sm,pm,xs,xt0,xt)
    nagents = par.nagents;

    crlb.rmt = zeros(1,par.nk);

    % --- INITIALIZE ---
    xr = xt0(1:2)-xs{1}; 
    H = sm.J(xr); H = H(1:2,1:2);
    P = blkdiag(inv(H'/sm.R*H),pm.v0^2*eye(2));
    for i = 2:nagents
        xr = xt0(1:2)-xs{i}; 
        H = sm.J(xr);
        P = inv(inv(P)+H'/sm.R*H);
    end
    
    % --- FILTER ---
    F = pm.F; Q = pm.Q;
    for k = par.kvec
        
        % TU:
        P = F*P*F' + Q;

        % MU:
        FI = inv(P);
        for i = 1:nagents
            xr = xt(1:2,k)-xs{i};
            H = sm.J(xr);
            FI = FI + H'/sm.R*H;
        end
        P = inv(FI);

        crlb.rmt(k) = sqrt(trace(P));
    end   
end

function stats = compute_nees_statistics(par,data,xt,iagent)
    idx = 1:4; m = length(idx);
    stats.M = par.M;
    stats.m = m;
    stats.idx = idx;

    s.Xi = zeros(m,m,par.nk);
    s.nees = zeros(1,par.nk);
    s.lambdamean = zeros(1,par.nk);
    s.lambdamin = zeros(1,par.nk);
    s.lambdamax = zeros(1,par.nk);

    % --- LKF ---
    stats.lkf = s;
    for k = par.kvec
        Xi = zeros(m);
        for imc = 1:par.M
            xtilde = data.lkf{imc,iagent}.xhat(idx,k) - xt(idx,k);
            L = chol(data.lkf{imc,iagent}.P(idx,idx,k),'lower'); 
            nu = L\xtilde;
            Xi = Xi + nu*nu';
        end
        Xi = Xi/par.M;
        lambda = eig(Xi);
        stats.lkf.Xi(:,:,k) = Xi;
        stats.lkf.nees(k) = trace(Xi);
        stats.lkf.lambdamean(k) = stats.lkf.nees(k)/m;
        stats.lkf.lambdamin(k) = min(lambda);
        stats.lkf.lambdamax(k) = max(lambda);
    end

    % --- NKF ---
    stats.nkf = s;
    for k = par.kvec
        Xi = zeros(m);
        for imc = 1:par.M
            xtilde = data.nkf{imc,iagent}.xhat(idx,k) - xt(idx,k);
            L = chol(data.nkf{imc,iagent}.P(idx,idx,k),'lower'); 
            nu = L\xtilde;
            Xi = Xi + nu*nu';
        end
        Xi = Xi/par.M;
        lambda = eig(Xi);
        stats.nkf.Xi(:,:,k) = Xi;
        stats.nkf.nees(k) = trace(Xi);
        stats.nkf.lambdamean(k) = stats.nkf.nees(k)/m;
        stats.nkf.lambdamin(k) = min(lambda);
        stats.nkf.lambdamax(k) = max(lambda);
    end

    % --- CI ---
    stats.ci = s;
    for k = par.kvec
        Xi = zeros(m);
        for imc = 1:par.M
            xtilde = data.ci{imc,iagent}.xhat(idx,k) - xt(idx,k);
            L = chol(data.ci{imc,iagent}.P(idx,idx,k),'lower'); 
            nu = L\xtilde;
            Xi = Xi + nu*nu';
        end
        Xi = Xi/par.M;
        lambda = eig(Xi);
        stats.ci.Xi(:,:,k) = Xi;
        stats.ci.nees(k) = trace(Xi);
        stats.ci.lambdamean(k) = stats.ci.nees(k)/m;
        stats.ci.lambdamin(k) = min(lambda);
        stats.ci.lambdamax(k) = max(lambda);
    end

    % --- LE ---
    stats.le = s;
    for k = par.kvec
        Xi = zeros(m);
        for imc = 1:par.M
            xtilde = data.le{imc,iagent}.xhat(idx,k) - xt(idx,k);
            L = chol(data.le{imc,iagent}.P(idx,idx,k),'lower'); 
            nu = L\xtilde;
            Xi = Xi + nu*nu';
        end
        Xi = Xi/par.M;
        lambda = eig(Xi);
        stats.le.Xi(:,:,k) = Xi;
        stats.le.nees(k) = trace(Xi);
        stats.le.lambdamean(k) = stats.le.nees(k)/m;
        stats.le.lambdamin(k) = min(lambda);
        stats.le.lambdamax(k) = max(lambda);
    end
end

function stats = compute_Wishart_statistics(m,n,p)
    alpha = 1-p;
    stats.lambdamin_ev = ev_lambdamin_Wishart(m,n)/n;
    stats.lambdamin_cdf = inv_cdf_lambdamin_Wishart(m,n,alpha/2)/n;
    stats.lambdamax_ev = ev_lambdamax_Wishart(m,n)/n;
    stats.lambdamax_cdf = inv_cdf_lambdamax_Wishart(m,n,(1-alpha/2))/n;
end

function plot_performance_stats(kidx,par,stats,crlb,clr)

    lw = 1;
    kvec = par.kvec(kidx);
    c = crlb.rmt(kidx);

    % --- RMT ---
    subplot(1,2,1);hold on
    hlkf = plot(kvec,stats.lkf.rmt(kidx)./c,'-'); hlkf.Color = clr.lkf; hlkf.LineWidth = lw;
    hnkf = plot(kvec,stats.nkf.rmt(kidx)./c,'-'); hnkf.Color = clr.naive; hnkf.LineWidth = lw;
    hci = plot(kvec,stats.ci.rmt(kidx)./c,'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hle = plot(kvec,stats.le.rmt(kidx)./c,'-'); hle.Color = clr.le; hle.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('RMT'); 
    box on

    % --- RMSE ---
    subplot(1,2,2);hold on
    hlkf = plot(kvec,stats.lkf.rmse(kidx)./c,'-'); hlkf.Color = clr.lkf; hlkf.LineWidth = lw;
    hnkf = plot(kvec,stats.nkf.rmse(kidx)./c,'-'); hnkf.Color = clr.naive; hnkf.LineWidth = lw;
    hci = plot(kvec,stats.ci.rmse(kidx)./c,'-'); hci.Color = clr.ci; hci.LineWidth = lw;
    hle = plot(kvec,stats.le.rmse(kidx)./c,'-'); hle.Color = clr.le; hle.LineWidth = lw;
    xlabel('$k$','interpreter','latex'); ylabel('RMSE');
    box on
end

function h = plot_nees_statistics(kidx,par,stats,clr)
    lw = 1; kvec = par.kvec(kidx);
    lmin = stats.lambdamin(kidx); lmax = stats.lambdamax(kidx); lambdamean = stats.lambdamean(kidx);
    h = plot(kvec,lambdamean,'-'); h.Color = clr; h.LineWidth = lw;
end

function h = plot_nees_region(kidx,par,stats,clr)
    kvec = par.kvec(kidx); lmin = stats.lambdamin(kidx); lmax = stats.lambdamax(kidx);
    opacity = 0.25;  
    xr = [kvec fliplr(kvec)];
    yr = [lmin fliplr(lmax)];    
    h = fill(xr,yr,clr);
    h.FaceAlpha = opacity;
    h.EdgeColor = clr;
    h.EdgeAlpha = opacity;
end

function h = plot_Wishart_statistics(kidx,par,stats,clr)
    lw = 1;
    kvec = par.kvec(kidx);
    lmax_ev = stats.lambdamax_ev*ones(size(kvec));
    lmax_cdf = stats.lambdamax_cdf*ones(size(kvec));
    h = plot(kvec,lmax_ev,'-'); h.Color = clr; h.LineWidth = lw;
    h = plot(kvec,lmax_cdf,'--'); h.Color = clr; h.LineWidth = lw;
end




% --- MISC ----------------------------------------------------------------
function idx = get_nonzero_indices(v)
    idx = zeros(1,length(v)); 
    for i = length(v):-1:1
        if v(i) == 0; idx(i) = []; 
        else; idx(i) = i;
        end
    end
    idx = num2cell(idx);
end

function idx = get_index_vector_diff(n,idx_remove)
    idx = 1:n;
    idx(idx_remove) = [];
end

function varargout = get_estimate_struct(varargin)
    s.xhat = []; s.P = [];
    if nargin > 0; s.xhat = varargin{1}; end
    if nargin > 1; s.P = varargin{2}; end
    for i = 1:nargout; varargout{i} = s; end
end

function s = get_record_struct(par)
    s.lkf = cell(par.M,par.nagents);
    s.nkf = s.lkf;
    s.ci = s.lkf;
    s.le = s.lkf;
    est.xhat = zeros(4,par.nk);
    est.P = zeros(4,4,par.nk);
    for imc = 1:par.M
        for i = 1:par.nagents
            s.lkf{imc,i} = est;
        end
    end
end

function s = get_dl_schedule(par)
    ntx = floor(par.nk/par.nagents); nrem = mod(par.nk,par.nagents); s = [];
    for i = 1:ntx; s = [s eye(par.nagents)]; end
    s = [s eye(par.nagents,nrem)];
end

function [x,varargout] = get_central_motion(ac,x0,N)
    p0 = x0(1:2); head0 = atan2(x0(4),x0(3));
    v = norm(x0(3:4));
    r = v^2/ac;
    L = N*v;
    ang_step = L/(r*N); 
    ang = 0;
    x = zeros(4,N);
    for k = 1:N
        ang = ang + ang_step;
        x(1:2,k) = r*[sin(ang) ; (1-cos(ang))];
        x(3:4,k) = v*[cos(ang) ; sin(ang)];
    end
    T = [cos(head0) -sin(head0) ; sin(head0) cos(head0)];
    x(1:2,:) = T*x(1:2,:) + p0;
    x(3:4,:) = T*x(3:4,:);
end

function add_plot_description(title_str)
    title(title_str)
    xlabel('$k$','interpreter','latex')
    ylabel('$\lambda$','interpreter','latex')
    box on
end