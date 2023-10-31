function params = tracker_parameters(varargin)
% --- tracker_parameters() ------------------------------------------------
% Returns parameters for a target tracker.
%
% 2023-10-30 Robin Forsling

% WITH DESCRIPTIONS
params.time = 'current simulated time';
params.all_tracks_init = 'set to true if all local tracks are initiated';
params.ntargets_in = 'number of targets input to the tracker';
params.ntracks_rx = 'number of received tracks';

if nargin > 0 && check_initiate_with_description(varargin{1}); return; 
else; params = set_struct_fields_empty(params);
end

% DEFAULTS
params.time = [];
params.all_tracks_init = 0;
params.ntargets_in = 0;
params.ntracks_rx = 0;