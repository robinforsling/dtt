function params = association_model_parameters(varargin)
% --- association_model_parameters() --------------------------------------
% Returns parameters for the association logic.
%
% 2023-10-30 Robin Forsling

% WITH DESCRIPTIONS
params.method = 'method used for association, e.g.: bypass, nn and gnn';
params.algorithm = 'algorithm used to solve assignment problem, e.g.: auction';
params.gate_prob = 'gating probability, where domain is [0,1), e.g.: 0.995';
params.llr = 'set to true if LLR should be used instead of Mahalanobis distance';

if nargin > 0 && check_initiate_with_description(varargin{1}); return; end

% DEFAULTS
params.method = def_association_method('bypass');
params.algorithm = [];
params.gate_prob = 0.995;
params.llr = 0;