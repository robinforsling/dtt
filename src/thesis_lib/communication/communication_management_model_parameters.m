function params = communication_management_model_parameters(varargin)
% --- communication_management_model_parameters() -------------------------
% Returns parameters for the communication management logic.
%
% 2023-10-30 Robin Forsling

% WITH DESCRIPTIONS
params.method = 'method used for selecting information, e.g.: pco, gevo, or dca';
params.ndirections = 'number of directions to select';
params.loss_function = 'loss function used for selection';

% DEFAULTS
params.method = def_communication_management_technique('gevo');
params.ndirections = 1;
params.loss_function = def_dimred_loss_function('ci'); 


