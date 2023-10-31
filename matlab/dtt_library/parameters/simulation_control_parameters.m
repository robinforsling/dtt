function cntrl = simulation_control_parameters(varargin)
% --- simulation_control_parameters() -------------------------------------
% Returns simulation control parameters.
%
% 2023-10-30 Robin Forsling

% Datalink
dl.req_estimates = 'datalink: request tracker estimates for transmission';

% Cheating
cheat.globally_known_est = 'cheat: global knowledge about covariances';
cheat.tracks = 'cheat: tracks to have global knowledge about';

% Set control struct
cntrl.dl = dl;
cntrl.cheat = cheat;
