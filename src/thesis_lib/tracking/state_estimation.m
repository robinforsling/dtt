function agent = state_estimation(time,agent,varargin)
% --- state_estimation() --------------------------------------------------
% Run measurement-to-track association and state estimation.
%
% VARARGIN
%   cntrl = varargin{1}: control inputs
%   targets = varargin{2}: input targets
%
% 2023-10-30 Robin Forsling

[cntrl,params,targets] = read_inputs(time,agent,nargin,varargin);

% GET
tracks = agent.tracks;
pm = agent.process_model; sm = agent.sensor_model; am = agent.association_model; fm = agent.fusion_model; dm = agent.datalink_model; 

% RUN TRACKER
meas = generate_local_measurements(sm,targets,params);
A = measurement_to_track_association(am,sm,tracks,meas,params);
tracks = kf_time_update(pm,tracks,params);
tracks = ekf_measurement_update(sm,tracks,meas,params,A);
%tracks = generate_datalink_estimates(cntrl,dm,fm,tracks,params);
[tracks,params.all_tracks_init] = track_initialization(pm,sm,tracks,meas,params);

% SET 
agent.params = params;
agent.tracks = tracks;

end
% -------------------------------------------------------------------------


% --- READ INPUTS ---------------------------------------------------------
function [cntrl,params,targets_in] = read_inputs(time,agent,NARG,VARGIN)

    % CONTROL
    cntrl = [];
    targets_in = [];
    
    % PARAMETERS
    params = agent.params;
    params.time = time;
    params.ntargets_in = 0;
    
    % HANDLE VARARGINS
    if NARG > 1; cntrl = VARGIN{1}; end 
    if NARG > 2; targets_in = VARGIN{2}; params.ntargets_in = length(targets_in); end
    
end
% -------------------------------------------------------------------------


