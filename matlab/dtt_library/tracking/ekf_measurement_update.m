function tracks = ekf_measurement_update(sensor_model,tracks,meas,params,varargin)
% --- ekf_measurement_update() --------------------------------------------
% Extended Kalman filter measurement update of tracks according using
% specified sensor model, tracks and input measurements. 
%
% VARARGIN
%   A = varargin{1}: assignment matrix (measurement to track assocation)
%
% 2023-10-30 Robin Forsling


% CHECK & HANDLE INPUTS
if isempty(meas); return; else; ntracks = length(tracks); nmeas = length(meas); end
if iscell(tracks); is_cell = 1; else; is_cell = 0; tracks = {tracks}; end
time = params.time;


% UPDATE TRACKS
m = sensor_model.m; h = sensor_model.h; J = sensor_model.J; R = sensor_model.R; xs = sensor_model.pos;
if ntracks > 1
    if nargin > 4; A = varargin{1}; else; error('multiple tracks require assignment matrix A as input'); end
    for i = 1:ntracks 
        if ~isnan(A(i,2)) && tracks{i}.init
            y = meas{A(i,2)}.y;
            tracks{i} = measurement_update(tracks{i},xs,h,J,R,y);
            tracks{i}.time = time;
        end
    end
else 
    if iscell(meas); y = meas{1}.y; else; y = meas.y; end
    tracks = tracks{1};
    if tracks.init
        tracks = measurement_update(tracks,xs,m,h,J,R,y);
        tracks.time.last_mu = time;
    end
    if is_cell; tracks = {tracks}; end
end

end


% --- MEASUREMENT UPDATE --------------------------------------------------
function track = measurement_update(track,xs,m,h,J,R,y)
    n = length(track.xhat);
    xrel = track.xhat(1:length(xs))-xs;
    H = J(xrel); H = [H zeros(m,n-size(H,2))];
    yhat = h(xrel); 
    ytilde = y-yhat; ytilde(2) = azimuth_modulus(ytilde(2));
    S = H*track.P*H' + R;
    K = track.P*H'/S;
    track.xhat = track.xhat + K*ytilde;
    track.P = (eye(n)-K*H)*track.P;
end


% --- AZIMUTH MODULUS -----------------------------------------------------
function a = azimuth_modulus(a)
    a = mod(a+pi,2*pi)-pi;
end

