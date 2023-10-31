function example_neglecting_correlations
% --- example_neglecting_correlations() -----------------------------------
% Example 3.9: Neglecting Correlations
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- TARGET --------------------------------------------------------------
target.pos = get_target_trajectory(2);


% --- SIMULATION PARAMETERS -----------------------------------------------
sim.M = 50;
sim.N = size(target.pos,2);
sim.process_model = 'cam';                                                  % cvm or cam


% --- AGENTS --------------------------------------------------------------
agent = get_agents(sim);

agent{1}.sensor.pos = [30;50];
agent{1}.sensor.R = [400 -30 ; -30 50];

agent{2}.sensor.pos = [200;30];
agent{2}.sensor.R = [400 30 ; 30 50];

for i = 1:2
    agent{i}.process.P0(1:2,1:2) = agent{i}.sensor.R;
end


% --- DATA RECORD ---------------------------------------------------------
data = cell(2,1); 
C = cell(sim.M,1); 
n = agent{1}.process.n;
for m = 1:sim.M; C{m}.xhat = zeros(n,sim.N); C{m}.P = zeros(n,n,sim.N); end
for i = 1:2; data{i}.lkf = C; data{i}.kf = C; data{i}.ci = C; end


% --- MAIN LOOP -----------------------------------------------------------
for m = 1:sim.M

    % --- INITIALIZE ----------------------------------------------------------
    for i = 1:2
        z = generate_measurement(target.pos(:,1),agent{i});
        x0 = [z ; zeros(agent{i}.process.n-2,1)]; 
        P0 = agent{i}.process.P0;
        agent{i}.est.lkf.xhat = x0; agent{i}.est.lkf.P = P0;
        agent{i}.est.kf.xhat = x0; agent{i}.est.kf.P = P0;
        agent{i}.est.ci.xhat = x0; agent{i}.est.ci.P = P0;
    end

    data = save_to_data_record(data,agent,m,1);

    for k = 2:sim.N
        
        % --- STATE ESTIMATION --------------------------------------------
        for i = 1:2
            z = generate_measurement(target.pos(:,k),agent{i});
            agent{i} = run_kalman_filter(z,agent{i});
        end
    
    
        % --- COMMUNICATION MANAGEMENT ------------------------------------
        if mod(k,2) == 0; itx = 1; irx = 2;
        else; itx = 2; irx = 1; 
        end
    
        est_tx = agent{itx}.est;
    
    
        % --- TRACK FUSION ------------------------------------------------
        agent{irx} = fuse_tracks(est_tx,agent{irx});
    
        data = save_to_data_record(data,agent,m,k);

    end
end

plot_estimates(target,data,sim)


end




% --- FUNCTIONS -----------------------------------------------------------
function z = generate_measurement(xt,agent)
    L = chol(agent.sensor.R,'lower');
    z = xt + L*randn(2,1);
end

function agent = run_kalman_filter(z,agent)

    F = agent.process.F; Q = agent.process.Q; 
    R = agent.sensor.R; H = agent.sensor.H;
    n = agent.process.n;

    % --- LKF -------------------------------------------------------------
    xhat = agent.est.lkf.xhat; P = agent.est.lkf.P;
    
    xhat = F*xhat;
    P = F*P*F' + Q;

    ztilde = z - H*xhat;
    S = H*P*H' + R;
    K = P*H'/S;
    agent.est.lkf.xhat = xhat + K*ztilde;
    agent.est.lkf.P = (eye(n)-K*H)*P;
    
    % --- KF --------------------------------------------------------------
    xhat = agent.est.kf.xhat; P = agent.est.kf.P;
    
    xhat = F*xhat;
    P = F*P*F' + Q;

    ztilde = z - H*xhat;
    S = H*P*H' + R;
    K = P*H'/S;
    agent.est.kf.xhat = xhat + K*ztilde;
    agent.est.kf.P = (eye(n)-K*H)*P;

    % --- CI --------------------------------------------------------------
    xhat = agent.est.ci.xhat; P = agent.est.ci.P;

    xhat = F*xhat;
    P = F*P*F' + Q;

    ztilde = z - H*xhat;
    S = H*P*H' + R;
    K = P*H'/S;
    agent.est.ci.xhat = xhat + K*ztilde;
    agent.est.ci.P = (eye(n)-K*H)*P;
end

function agent = fuse_tracks(est_tx,agent)
    agent.est.kf = kalman_fuser(est_tx.kf,agent.est.kf);
    agent.est.ci = covariance_intersection(est_tx.ci,agent.est.ci);
end

function est = kalman_fuser(est1,est2)
    y1 = est1.xhat; R1 = est1.P; I1 = inv(R1);
    y2 = est2.xhat; R2 = est2.P; I2 = inv(R2);
    
    est.P = inv(I1 + I2);
    est.xhat = est.P*(I1*y1 + I2*y2);
end

function est = covariance_intersection(est1,est2)
    y1 = est1.xhat; R1 = est1.P; I1 = inv(R1);
    y2 = est2.xhat; R2 = est2.P; I2 = inv(R2);

    f = @(w) trace(inv(w*I1 + (1-w)*I2));
    w = fminbnd(f,0,1);

    est.P = inv(w*I1 + (1-w)*I2);
    est.xhat = est.P*(w*I1*y1 + (1-w)*I2*y2);
end

function data = save_to_data_record(data,agent,m,k)
    for i = 1:2
        data{i}.lkf{m}.xhat(:,k) = agent{i}.est.lkf.xhat;
        data{i}.lkf{m}.P(:,:,k) = agent{i}.est.lkf.P;
        data{i}.kf{m}.xhat(:,k) = agent{i}.est.kf.xhat;
        data{i}.kf{m}.P(:,:,k) = agent{i}.est.kf.P;
        data{i}.ci{m}.xhat(:,k) = agent{i}.est.ci.xhat;
        data{i}.ci{m}.P(:,:,k) = agent{i}.est.ci.P;        
    end
end 




% --- PLOT FUNCTIONS ------------------------------------------------------
function plot_estimates(target,data,sim)
    ia = 1;
    clr = get_thesis_colors;
    lwe = 0.5; lwt = 2;

    % --- FIGURE ----------------------------------------------------------
    figure(1);clf;hold on
    for m = 1:sim.M 
        hlkf = data{ia}.lkf{m}.xhat; xkf = data{ia}.kf{m}.xhat; xci = data{ia}.ci{m}.xhat;
        hkf = plot(xkf(1,:),xkf(2,:),'k'); hkf.Color = clr.kf; hkf.LineWidth = lwe;
        hci = plot(xci(1,:),xci(2,:),'k'); hci.Color = clr.ci; hci.LineWidth = lwe;
    end
    ht = plot(target.pos(1,:),target.pos(2,:),'k-');  ht.Color = clr.black; ht.LineWidth = lwt;
    ht0 = plot(target.pos(1,1),target.pos(2,1),'ko'); ht0.Color = clr.black; ht0.MarkerFaceColor = clr.black;

    box on; axis equal
    legend(gca,[ht ht0 hkf hci],'tagret trajectory','target initial position','KF','CI','location','southwest')
    
    set_fontsize_all(14)
end



% --- OTHER FUNCTIONS -----------------------------------------------------
function xt = get_target_trajectory(ID)
    x0 = [0;0];
    switch ID
        case 1
            a = straight_segment_fixed_length(x0,0,240,13);
            b = circle_arc_fixed_length_right(a,63,100,6);
            c = straight_segment_fixed_length(b,60,4);
            merged = merge_motion_primitives(a,b,c);
    
        case 2
            a = straight_segment_fixed_length(x0,0,140,8);
            b = circle_arc_fixed_length_right(a,60,100,6);
            c = circle_arc_fixed_length_left(b,60,100,6);
            d = straight_segment_fixed_length(c,100,6);
            merged = merge_motion_primitives(a,b,c,d);

        case 3
            a = straight_segment_fixed_length(x0,0,240,13);
            b = circle_arc_fixed_length_right(a,50,160,9);
            c = straight_segment_fixed_length(b,100,6);
            merged = merge_motion_primitives(a,b,c);
    end
    xt = merged.X;
end

function agent = get_agents(sim)
    a.est = get_estimate_struct;
    a.process = [];
    a.sensor.pos = []; a.sensor.H = []; a.sensor.R = [];
    agent = cell(2,1);
    for i = 1:2
        agent{i} = a;
        if strcmpi(sim.process_model,'cvm')
            agent{i}.process = get_cvm;
            agent{i}.sensor.H = [eye(2) zeros(2)];
        elseif strcmpi(sim.process_model,'cam')
            agent{i}.process = get_cam;
            agent{i}.sensor.H = [eye(2) zeros(2,4)];
        else
            error('unknown process model...')
        end
    end
end

function est = get_estimate_struct()
    est.kf.xhat = [];
    est.kf.P = [];
    est.ci.xhat = [];
    est.ci.P = [];
end

function pm = get_cvm()
    T = 1;
    q = 4;
    pm.n = 4; pm.T = T; pm.q = q;
    pm.F = [eye(2) T*eye(2) ; zeros(2) eye(2)];
    pm.Q = q^2*[T^3*eye(2)/3 T^2*eye(2)/2 ; T^2*eye(2)/2 T^2*eye(2)];
    pm.vmax = 30;
    pm.P0 = [];
end

function pm = get_cam()
    T = 1;
    q = 1.0;
    pm.n = 6; pm.T = T; pm.q = q;
    pm.F = [eye(2) T*eye(2) T^2*eye(2)/2 ; zeros(2) eye(2) T*eye(2) ; zeros(2,4) eye(2)];
    pm.Q = q^2*[T^5*eye(2)/20 T^4*eye(2)/8 T^3*eye(2)/6 ; ...
                T^4*eye(2)/8 T^3*eye(2)/3 T^2*eye(2)/2 ; ...
                T^3*eye(2)/6 T^2*eye(2)/2 T*eye(2)];
    pm.vmax = 30;
    pm.amax = 2;
    pm.P0 = blkdiag(zeros(2),pm.vmax^2*eye(2),pm.amax^2*eye(2));
end
