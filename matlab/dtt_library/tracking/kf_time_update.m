function tracks = kf_time_update(process_model,tracks,params)
% --- kf_time_update() ----------------------------------------------------
% Kalman filter prediction of track according to linear process model given
% by process_model.
%
% 2023-10-30 Robin Forsling


% HANDLE INPUTS
if iscell(tracks); is_cell = 1; else; is_cell = 0; tracks = {tracks}; end
time = params.time;

if isempty(process_model.F) 
    Ts = model.Ts; q = model.q;
    F = process_model.Fbasis(Ts); Q = q^2*process_model.Qbasis(Ts);
else
    F = process_model.F; Q = process_model.Q;
end

% RUN PREDICTION
for i = 1:length(tracks)
    if tracks{i}.init
        tracks{i}.xhat = F*tracks{i}.xhat;
        tracks{i}.P = F*tracks{i}.P*F' + Q;
        tracks{i}.time.last_tu = time;
        tracks{i}.common_est.xhat = F*tracks{i}.common_est.xhat;
        tracks{i}.common_est.P = F*tracks{i}.common_est.P*F' + Q;
    end
end

% CONVERT BACK TO STRUCT
if ~is_cell; tracks = tracks{1}; end

