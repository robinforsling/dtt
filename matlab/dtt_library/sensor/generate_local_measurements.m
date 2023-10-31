function meas = generate_local_measurements(sm,targets,params)
% --- generate_local_measurements() ---------------------------------------
% Generate measurement array for a multitarget tracker.
%
% 2023-10-30 Robin Forsling

% HANDLE INPUTS
if isempty(targets); meas = []; return; else; ntargets = length(targets); meas = cell(ntargets,1); end
time = params.time;

m = sm.m; h = sm.h; R = sm.R; L = chol(R,'lower');
PD = sm.PD; xs = sm.pos; ori = sm.ori; fov = sm.fov; r_lim = sm.r_lim;

% CONSTRUCT MEASUREMENTS
for i = 1:ntargets
    xt = targets{i}.pos; xrel = xt-xs;
    meas{i} = create_measurement_struct;
    if check_detection(PD,ori,fov,xrel)
        y = h(xrel) + L*randn(m,1); 
        y(1) = max([y(1) r_lim(1)]);
        y(2) = azimuth_modulus(y(2));
        meas{i}.y = y;
        meas{i}.ID = targets{i}.ID;
        meas{i}.detection = 1;
        meas{i}.time = time;
    else
        meas{i}.detection = 0;
    end
end

end
% -------------------------------------------------------------------------


% --- CHECK DETECTION -----------------------------------------------------
function target_detected = check_detection(PD,ori,fov,xrel)
    target_detected = 1;
    b = atan2(xrel(2),xrel(1)); b = mod(b+pi,2*pi); ori = mod(ori+pi,2*pi);
    if ~is_target_in_fov(b,ori,fov); target_detected = 0; return; end
    if rand > PD; target_detected = 0; end
end
% -------------------------------------------------------------------------


% --- CHECK IF TARGET IS IN FOV -------------------------------------------
function in_fov = is_target_in_fov(b,ori,fov)  
    in_fov = 0; if abs(b-ori) <= fov/2; in_fov = 1; end
end
% -------------------------------------------------------------------------


% --- AZIMUTH MODULUS -----------------------------------------------------
function a = azimuth_modulus(a)
    a = mod(a+pi,2*pi)-pi;
end
% -------------------------------------------------------------------------
