function out_data = simenv(scen_par)
% --- simenv() ------------------------------------------------------------
% Simulation environment for decentralized target tracking (DTT).
%
% 2024-01-18 Robin Forsling


% --- SETUP SIM -----------------------------------------------------------

% PARAMETERS
scen = load_scen(scen_par);

nagents = scen.nagents; ntargets = scen.ntargets;
nx = scen.agents{1}.process_model.nx;
ncoord = scen.par.ncoord;
M = scen.par.M;
N = scen.par.N;

% DATA REC
data = cell(nagents,1);
for i = 1:nagents
    data{i}.xhat = cell(M,1);
    data{i}.P = cell(M,1);
end


% --- MISC ----------------------------------------------------------------
k0 = 1;
k25 = round(0.25*M);
k50 = round(0.5*M);
k75 = round(0.75*M);
k100 = M;


% --- MAIN LOOP -----------------------------------------------------------
for imc = 1:M

    switch imc
        case k0;  fprintf('    %s: 0%%\n',get_datetime)
        case k25; fprintf('    %s: 25%%\n',get_datetime)
        case k50; fprintf('    %s: 50%%\n',get_datetime)  
        case k75; fprintf('    %s: 75%%\n',get_datetime)   
    end

    agents = scen.agents;
    targets = scen.targets;

    XHAT = cell(nagents,1);
    PP = cell(nagents,1);
    for i = 1:nagents
        XHAT{i} = zeros(nx,N);
        PP{i} = zeros(nx,nx,N);
    end

    for k = 1:N
    
        % Update targets
        for itgt = 1:ntargets
            targets{itgt}.pos = scen.targets{itgt}.pos(:,k);
            targets{itgt}.vel = scen.targets{itgt}.vel(:,k);
            targets{itgt}.acc = scen.targets{itgt}.acc(:,k);
        end

        % Communication/network dynamics
        dl_schedule = scen.dl.schedule;
        if ndims(scen.dl.topology) == 3; dl_topology = scen.dl.topology(:,:,k);
        else; dl_topology = scen.dl.topology;
        end
        if sum(diag(dl_topology)) > 0; error('error in dl topology...'); end
        %dl_schedule = sum(dl_topology,2); dl_schedule(dl_schedule > 0) = 1;

        % Simulation control parameters
        cntrl = simulation_control_parameters;
        cntrl.cheat.globally_known_est = scen.cntrl.globally_known_est;

        % Run track filter
        for i = 1:nagents
            agents{i} = state_estimation(k,agents{i},cntrl,targets);
        end
    
        % Run agent subsystems
        for i = 1:nagents
            if scen.cntrl.globally_known_est && dl_schedule(i,k) % Realization of "globally known Ri cheat"
                if sum(dl_topology(i,:)) > 1
                    error('this case is not implemented')
                else
                    for j = 1:nagents
                        if dl_topology(i,j) && i ~= j; cntrl.cheat.tracks = agents{j}.tracks; end
                    end
                end
            end
            cntrl.dl.req_estimates = dl_schedule(i,k);
            agents{i} = communication_management(k,agents{i},cntrl,targets);
        end
    
        % Simulate communication
        dl_xmit = cell(nagents,1);
        n_xmit = 0;
        for i = 1:nagents
            if dl_schedule(i,k) 
                dl_xmit{i}.tracks = agents{i}.tracks; 
                n_xmit = n_xmit + 1;
            else 
                dl_xmit{i}.tracks = [];
            end
        end
    
        % Run track fusers
        for i = 1:nagents
            cntrl = []; dl_tracks = cell(n_xmit,1);
            irpt = 1;
            for j = 1:nagents
                if dl_schedule(j,k) && dl_topology(j,i) && i ~= j
                    dl_tracks{irpt} = dl_xmit{j}.tracks;
                    irpt = irpt + 1;
                end
            end
            agents{i} = track_fusion(k,agents{i},cntrl,dl_tracks);
        end

        %  Save data
        for i = 1:nagents
            XHAT{i}(:,k) = agents{i}.tracks{1}.xhat;
            PP{i}(:,:,k) = agents{i}.tracks{1}.P;
        end
    end

    for i = 1:nagents
        data{i}.xhat{imc} = XHAT{i};
        data{i}.P{imc} = PP{i};
    end

    if imc == k100; fprintf('    %s: 100%%\n',get_datetime); end
end

out_data.scen_par = scen.par;
out_data.sensor_pos = scen.sensor_pos;
out_data.dl = scen.dl;
out_data.cntrl = scen.cntrl;
out_data.nagents = nagents;
out_data.ntargets = ntargets;
out_data.agents = agents;
out_data.targets = scen.targets;
out_data.rec = data;

