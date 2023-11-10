function s = default_simulation_parameters
% --- default_simulation_parameters() -------------------------------------
% Default simulation parameters.
%
% 2023-11-07

% Basic stuff
s.ID = 0;
s.m = 1;
s.comm_mgmt_method = 4;
s.dimred_loss = 4;
s.fus_method = 4;
s.M = 100;
s.ncoord = 2;
s.dl_T = 2;
s.sigma_r = 1000;
s.sigma_az = 0.1*d2r;

% Control
c.globally_known_est = 0;

s.cntrl = c;