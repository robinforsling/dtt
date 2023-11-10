function track = track_parameters(varargin)
% --- track_parameters() --------------------------------------------------
% Returns struct of track parameters.
%
% 2023-10-30 Robin Forsling

time.init = [];                                                             % Time of initiation
time.last_tu = [];                                                          % Time of last time update
time.last_mu = [];                                                          % Time of last measurement update
time.last_fuse = [];                                                        % Time of last fusion

% MAIN ESTIMATE
track.ID = [];                                                              % Track ID
track.init = 0;                                                             % True if track has been initiated, otherwise false
track.xhat = [];                                                            % State estimate
track.P = [];                                                               % State estimate covariance
%track.H = [];
track.time = time;

% COMMON ESTIMATE
%track.common_est.valid = 0;
track.common_est.xhat = [];
track.common_est.P = [];
%track.common_est.H = [];
%track.common_est.time = time;

track.dl_est.valid = 0;
track.dl_est.xhat = [];
track.dl_est.P = [];
track.dl_est.H = [];
track.dl_est.time = [];
track.dl_est.time_last_xmit = [];




