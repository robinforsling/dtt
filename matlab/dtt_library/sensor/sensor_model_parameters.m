function model = sensor_model_parameters(varargin)
% --- sensor_model_parameters() -------------------------------------------
% Creates parameters for a sensor model.
%
% 2023-10-30 Robin Forsling

init_with_descr = 0;
if nargin > 0; init_with_descr = check_initiate_with_description(varargin{1}); end

% With descriptions
model.sensor_type = 'type of sensor: polar or spherical';
model.pos = 'sensor position';
model.ori = 'sensor orientation in horizontal plane';
model.r_lim = 'sensor range limits';
model.fov = 'sensor field-of-view in horizontal plane';
model.m = 'measurement dimensionality: 2 or 3';          
model.h = 'measurement function';
model.J = 'Jacobian, derivative of h w.r.t. x';
model.R = 'measurement covariance';
model.PD = 'probability of detection';
model.beta_FA = 'false alarm rate';

% Set all fields empty
if ~init_with_descr; model = set_struct_fields_empty(model); end