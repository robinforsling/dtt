function [tracks,varargout] = track_initialization(process_model,sensor_model,tracks,meas,params)
% --- track_initialization() ----------------------------------------------
% Initialized tracks not already initiated.
%
% VARARGOUT
%   varargout{1} = all_tracks_init: true if all tracks are already initiated 
%
% 2023-10-30 Robin Forsling

% HANDLE INPUT
if params.all_tracks_init; if nargout > 1; varargout{1} = 1; end; return; end             % Return if all tracks already initiated
if ~iscell(tracks); tracks = {tracks}; end
if ~iscell(meas); meas = {meas}; end
time = params.time;

ntracks = length(tracks);
not_init_tracks = [];

for i = 1:ntracks
    if ~tracks{i}.init
        if meas{i}.detection
            tracks{i} = init_track(process_model,sensor_model,tracks{i},meas{i});
            tracks{i}.init = 1;
            tracks{i}.time.init = time;
        else
            not_init_tracks = [not_init_tracks i];
        end
    end
end

if nargout > 1; varargout{1} = isempty(not_init_tracks); end

end
% -------------------------------------------------------------------------


% --- INIT TRACK ----------------------------------------------------------
function track = init_track(pm,sm,track,meas)
    smax2 = max(sm.R(1,1),sm.R(2,2)*meas.y(1)^2); 
    I = eye(pm.ncoord);
    switch sm.m
        case 2; x = convert_pol_to_cart(meas.y);
        case 3; x = convert_sph_to_cart(meas.y);
        otherwise; error('initialize_tracks: invalid dimensionality')
    end
    if strcmpi(pm.model_type,'cv'); P0 = blkdiag(smax2*I,pm.v_max^2*I);
    elseif strcmpi(pm.model_type,'ca'); P0 = blkdiag(smax2*I,pm.v_max^2*I,pm.a_max^2*I);
    end
    track.xhat = [x+sm.pos;zeros(pm.nx-sm.m,1)];
    track.P = P0;
    track.common_est.xhat = track.xhat;
    track.common_est.P = 1.001*P0;
end
% -------------------------------------------------------------------------

