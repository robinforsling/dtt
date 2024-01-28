function s = default_simulation_parameters
% --- default_simulation_parameters() -------------------------------------
% Default simulation parameters.
%
% 2024-01-15

% Basic stuff:
s.ID = 0;
s.M = 100;
s.ncoord = 2;

% Sensor parameters:
s.sigma_r = 1000;
s.sigma_az = 0.1*d2r;

% Filter parameters:
s.q = 1;
s.Ts = 1;

% Track fusion parameters:
s.fus_method = 4;

% Datalink parameters:
s.m = 1;
s.dl_T = 2;
s.comm_mgmt_method = 4;
s.dimred_loss = 4;

% Simulation control:
c.globally_known_est = 0;
s.cntrl = c;