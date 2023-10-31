function agent = communication_management(time,agent,varargin)
% --- communication_management() ----------------------------------------------
% Run communication management subsystem.
%
% VARARGIN
%   cntrl = varargin{1}: control inputs
%   targets = varargin{2}: input targets
%
% 2023-10-30 Robin Forsling

[cntrl,params,~] = read_inputs(time,agent,nargin,varargin);

% GET
tracks = agent.tracks;
fm = agent.fusion_model; dm = agent.datalink_model;                         %pm = agent.process_model; sm = agent.sensor_model; am = agent.association_model;  

% RUN SUBSYSTEMS
tracks = generate_datalink_tracks(cntrl,dm,fm,tracks,params);

% SET
agent.tracks = tracks;

end


% --- READ INPUTS ---------------------------------------------------------
function [cntrl,params,targets_in] = read_inputs(time,tracker,NARG,VARGIN)

    % CONTROL
    cntrl = [];
    targets_in = [];
    
    % PARAMETERS
    params = tracker.params;
    params.time = time;
    params.ntargets_in = 0;
    
    % HANDLE VARARGINS
    if NARG > 1; cntrl = VARGIN{1}; end 
    if NARG > 2; targets_in = VARGIN{2}; params.ntargets_in = length(targets_in); end
    
end