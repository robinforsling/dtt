function track = datalink_track_parameters()
% --- datalink_track_parameters() -----------------------------------------
% Returns struct of datalink track parameters.
%
% 2023-10-30 Robin Forsling

track.ID = [];                                                              % Track ID
track.init = 0;                                                             % True if track has been initiated, otherwise false

track.xhat = [];                                                            % State estimate
track.P = [];                                                               % State estimate covariance
track.H = [];
