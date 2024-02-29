function scen = load_scen(scen_par)
% --- load_scen() ---------------------------------------------------------
% Loads requested scenario. Additional scenarios are specified below.
%
% 2024-01-18 Robin Forsling


% --- DEFAULTS ------------------------------------------------------------
Ts = 1;
q = 1;
nagents = 2;
N = [];

switch scen_par.ID

    % --- SCENARIO 1: 2 AGENTS --------------------------------------------
    case 1
        N = 18; q = 2;

        tgt_traj = circle_arc_fixed_length_left([3000;8000],-pi/6,5000,4000,N);
        
        XS = cell(nagents,1); 
        XS{1} = [-2000;1000]; XS{2} = [5000;0]; 
        agents = cell(nagents,1);
        pm.Ts = 1; pm.model_type = 'ca';  pm.q = q;
        pm = process_model_linear(pm);
        pm.v_max = 170; pm.a_max = 5; 

        dl.topology = [0 1 ; 1 0];
        dl.schedule = basic_communication_schedule(nagents,N);

        for i = 1:nagents
            sm = sensor_model_spherical([],scen_par.sigma_r,scen_par.sigma_az); sm.pos = XS{i};
            agents{i} = create_single_target_tracker(pm,sm);
        end


    % --- SCENARIO 2: 2 AGENTS --------------------------------------------
    case 2
        N = 18; q = 4;

        tgt_traj = circle_arc_fixed_length_left([3000;8000],-pi/6,5000,4000,N);
        
        XS = cell(nagents,1); 
        XS{1} = [-2000;1000]; XS{2} = [5000;0]; 
        agents = cell(nagents,1);
        pm.Ts = 1; pm.model_type = 'ca';  pm.q = q;
        pm = process_model_linear(pm);
        pm.v_max = 200; pm.a_max = 5; 

        dl.topology = [0 1 ; 1 0];
        dl.schedule = basic_communication_schedule(nagents,N);

        for i = 1:nagents
            sm = sensor_model_spherical([],scen_par.sigma_r,scen_par.sigma_az); sm.pos = XS{i};
            agents{i} = create_single_target_tracker(pm,sm);
        end


    % --- SCENARIO 3 ------------------------------------------------------
    case 3
        N = 12; q = 1;

        tgt_traj = circle_arc_fixed_length_left([3000;7000],-pi/6,10000,4000,N);
        
        XS = cell(nagents,1); 
        XS{1} = [-4000;8000]; XS{2} = [2500;-1000]; 
        agents = cell(nagents,1);
        pm.Ts = 1; pm.model_type = 'ca';  pm.q = q;
        pm = process_model_linear(pm);
        pm.v_max = 200; pm.a_max = 10; 

        dl.topology = [0 1 ; 1 0];
        dl.schedule = basic_communication_schedule(nagents,N);

        for i = 1:nagents
            sm = sensor_model_spherical([],scen_par.sigma_r,scen_par.sigma_az); sm.pos = XS{i};
            agents{i} = create_single_target_tracker(pm,sm);
        end


    % --- SCENARIO 4 ------------------------------------------------------
    case 4
        N = 20; nagents = 3; Ts = 1; q = 5;

        tgt_traj = circle_arc_fixed_length_right([8000;8000],pi/6,3000,3000,N);
        
        XS = cell(nagents,1); 
        XS{1} = [0;6000]; XS{2} = [1000;0]; XS{3} = [0;10000];
        agents = cell(nagents,1);
        pm.Ts = 1; pm.model_type = 'ca';  pm.q = 1;
        pm = process_model_linear(pm);
        pm.v_max = 100; pm.a_max = 10; 

        dl.topology = [0 1 0  ; 1 0 0 ; 0 1 0];
        dl.schedule = [basic_communication_schedule(2,N);zeros(1,N)]; dl.schedule(3,11:2:N) = 1;

        sigma = [1000 3*d2r ; 1000 3*d2r ; 1000 0.75*d2r];

        for i = 1:nagents
            sm = sensor_model_spherical([],sigma(i,1),sigma(i,2)); sm.pos = XS{i};
            agents{i} = create_single_target_tracker(pm,sm);
        end


    % --- SCENARIO 5: 4 AGENTS --------------------------------------------
    case 5
        N = 40; nagents = 4; Ts = 1; q = 5;

        tgt_traj = circle_arc_fixed_length_left([3000;5000],-pi/12,6000,4000,N);
        
        XS = cell(nagents,1); 
        XS{1} = [-1000;10000]; XS{2} = [0;0]; XS{3} = [15000;1000]; XS{4} = [5000;15000];
        agents = cell(nagents,1);
        pm = process_model_linear; pm.model_type = 'ca'; pm = process_model_linear(pm);
        pm.v_max = 200; pm.a_max = 10;

        dl.topology = [0 1 0 0 ; 0 0 1 0 ; 1 0 0 1 ; 1 0 0 0];
        dl.schedule = basic_communication_schedule(nagents,N);
        %dl.schedule = zeros(ns,N); for i = 1:ns; dl.schedule(1:4:N+(i-1)) = 1; end

        for i = 1:nagents
            sm = sensor_model_spherical([],scen_par.sigma_r,scen_par.sigma_az); sm.pos = XS{i};
            agents{i} = create_single_target_tracker(pm,sm);
        end


    % --- SCENARIO 10: NEES MATRIX APPLICATION ----------------------------   
    case 10
        nagents = 2; agents = cell(nagents,1);
        Ts = 1; 

        XS = cell(nagents,1); XS{1} = [-2000;1000]; XS{2} = [5000;0]; 
        a = straight_segment_fixed_length([1000;6000],-2*pi/5,200,3,Ts); 
        b = circle_arc_fixed_length_left(a,2000,1000,11,Ts);
        c = straight_segment_fixed_length(b,200,3,Ts);
        d = circle_arc_fixed_length_right(c,6000,1000,11,Ts);
        tgt_traj = merge_motion_primitives(a,b,c,d);
        %tgt_traj = straight_segment_fixed_length([3000;5000],-pi/12,4000,40,Ts);

        pm.Ts = Ts; pm.model_type = 'cv';  pm.q = scen_par.q;
        pm = process_model_linear(pm);
        pm.v_max = 200; pm.a_max = 5; 
        sm = sensor_model_spherical([],scen_par.sigma_r,scen_par.sigma_az);

        dl.topology = [0 1 ; 1 0];
        dl.schedule = zeros(nagents,N); %basic_communication_schedule(nagents,N);

        for i = 1:nagents
            sm.pos = XS{i};
            agents{i} = create_single_target_tracker(pm,sm);
        end


    % --- SCENARIO 11: NEES MATRIX APPLICATION ---------------------------- 
    % Filter tuning and track fusion design
    case 11
        nagents = 2; agents = cell(nagents,1);
        Ts = 1; 

        XS = cell(nagents,1); XS{1} = [-3000;1000]; XS{2} = [5000;0]; 
        tgt_traj = circle_arc_fixed_length_left([2000;6000],-2*pi/3,5000,3000,31,Ts); 

        pm.Ts = Ts; pm.model_type = 'cv';  pm.q = scen_par.q;
        pm = process_model_linear(pm);
        pm.v_max = 100; pm.a_max = 3; 
        sm = sensor_model_spherical([],scen_par.sigma_r,scen_par.sigma_az);

        dl.topology = [0 1 ; 1 0];
        dl.schedule = basic_communication_schedule(nagents,N);

        for i = 1:nagents
            sm.pos = XS{i};
            agents{i} = create_single_target_tracker(pm,sm);
        end

    % --- SCENARIO 12: NIS MATRIX APPLICATION -----------------------------
    % Target maneuver detection
    case 12 
        nagents = 1; agents = cell(nagents,1);
        Ts = 1;

        XS = cell(nagents,1); XS{1} = [0;0];
        pm.Ts = Ts; pm.model_type = 'cv'; pm.q = scen_par.q;
        pm.v_max = 100; pm.a_max = 3; 
        sm = sensor_model_spherical([],scen_par.sigma_r,scen_par.sigma_az);

        a = straight_segment_fixed_length([10000;10000],5*pi/4,1000,11,Ts); 
        b = circle_arc_fixed_length_left(a,500,1000,11,Ts);
        tgt_traj = merge_motion_primitives(a,b);



    otherwise 
        error('load_scen: unknown scenario ID')

end

if ~iscell(tgt_traj); tgt_traj = {tgt_traj}; end
ntargets = length(tgt_traj);
targets = cell(ntargets,1);
for itgt = 1:ntargets
    targets{itgt} = target_parameters;
    targets{itgt}.ID = 1;
    targets{itgt}.X = tgt_traj{itgt}.X;
    targets{itgt}.pos = tgt_traj{itgt}.X(1:2,:);
    targets{itgt}.vel = tgt_traj{itgt}.X(3:4,:);
    targets{itgt}.acc = tgt_traj{itgt}.X(5:6,:);
    targets{itgt}.ori = tgt_traj{itgt}.head;
end

if isempty(N); N = tgt_traj{1}.N; end


% --- PACK DATA -----------------------------------------------------------
scen.par.ID = scen_par.ID;
scen.par.M = scen_par.M;
scen.par.N = N;
scen.par.Ts = Ts;
scen.par.ncoord = scen_par.ncoord;
scen.sensor_pos = XS;
scen.ntargets = ntargets;
scen.nagents = length(agents);
scen.targets = targets;
scen.agents = agents;
scen.dl = dl;
scen.cntrl = scen_par.cntrl;

for i = 1:scen.ntargets
    scen.targets{i}.ID = i;
    for j = 1:scen.nagents; scen.agents{j}.tracks{i}.ID = i; end
end

for i = 1:scen.nagents
    scen.agents{i}.datalink_model.comm_mgmt.method = scen_par.comm_mgmt_method;
    scen.agents{i}.fusion_model.method = scen_par.fus_method;
    scen.agents{i}.datalink_model.comm_mgmt.ndirections = scen_par.m;
    scen.agents{i}.datalink_model.comm_mgmt.loss_function = scen_par.dimred_loss;
end

