function agent = track_fusion(time,agent,varargin)
% --- track_fusion() ---------------------------------------------------
% Run multitarget tracking fuser.
%
% VARARGIN
%   cntrl = varargin{1}: control inputs
%   tracks_rx = varargin{2}: cell array of tracks received
%   time = varargin{3}: simulated time (current time)
%
% 2023-10-30 Robin Forsling

[cntrl,params,tracks_rx] = read_inputs(time,agent,nargin,varargin);

% GET
tracks_local = agent.tracks;
fm = agent.fusion_model; am = agent.association_model;

% RUN FUSER
for irx = 1:length(tracks_rx)
    A = track_to_track_association(am,fm,tracks_local,tracks_rx{irx},params);
    tracks_local = fuse_tracks(fm,tracks_local,tracks_rx{irx},params,A);
end

% SET 
agent.params = params;
agent.tracks = tracks_local;

end
% -------------------------------------------------------------------------


% --- READ INPUTS ---------------------------------------------------------
function [cntrl,params,tracks_rx] = read_inputs(time,agent,NARG,VARGIN)

    % CONTROL
    cntrl = [];
    
    % PARAMETERS
    params = agent.params;
    params.time = time;
    params.ntracks_rx = 0;
    
    % HANDLE VARARGINS
    tracks_rx = [];
    if NARG > 1; cntrl = VARGIN{1}; end 
    if NARG > 2; tracks_rx = VARGIN{2}; params.ntracks_rx = length(tracks_rx); end
    
end
% -------------------------------------------------------------------------


