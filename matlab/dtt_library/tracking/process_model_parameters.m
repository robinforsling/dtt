function model = process_model_parameters()
% --- process_model_parameters() ------------------------------------------
% Creates parameters for a process model.
%
% 2023-10-30 Robin Forsling

% WITH DESCRIPTIIONS
model.model_type = 'type of process model: cp, cv or ca';           
model.ncoord = 'number of spatial coordinates, e.g.: 2 or 3';           
model.nx = 'state dimensionality, e.g.: 4, 6 or 9';             
model.Ts = 'sampling time, e.g.: 1s or set to empty to get F and Q as function handles';   
model.F = 'process model';
model.Q = 'process noise covariance';
model.q = 'process noise magnitude/standard deviation';
model.v_max = 'assumed maximum velocity of target (in each coordinate)';
model.a_max = 'assumed maximum acceleration of target (in each coordinate)';