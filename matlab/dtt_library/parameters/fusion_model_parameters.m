function params = fusion_model_parameters(varargin)
% --- fusion_model_parameters() -------------------------------------------
% Returns parameter struct for track fusion functionality.
%
% 2023-10-30 Robin Forsling

params.method = 'fusion method to use, e.g.: CI or naive';

% DEFAULT
params.method = def_fusion_method('ci');