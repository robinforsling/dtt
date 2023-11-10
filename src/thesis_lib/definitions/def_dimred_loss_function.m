function out_var = def_dimred_loss_function(varargin)
% --- def_dimred_loss_function(varargin) ----------------------------------
% Definition of loss function used in the communication management
% dimension-reduction.
%
% 2023-10-30 Robin Forsling

defs = {1, 'kf'; ...                                                        % KF
        2, 'kf-rce'; ...                                                    % KF-remove-common-estimate
        3, 'bsc'; ...                                                       % BSC
        4, 'ci'; ...                                                        % CI
        5, 'ici'; ...                                                       % ICI (not implemented yet)
        6, 'le'};                                                           % LE              

% HANDLE INPUTS
input_is_string = 0; input_is_scalar = 0;
if nargin > 0
    var = varargin{1};
    if isstring(var) || ischar(var); input_is_string = 1; 
    elseif isnumeric(var) && isscalar(var); input_is_scalar = 1; end
else
    error('too few inputs')
end

% CONVERT
ndef = size(defs,1); out_var = [];
if input_is_string
    for i = 1:ndef; if strcmpi(defs{i,2},var); out_var = defs{i,1}; end; end
    if isempty(out_var); error('unknown input string'); end
elseif input_is_scalar
    for i = 1:ndef; if defs{i,1} == var; out_var = defs{i,2}; end; end
    if isempty(out_var); error('unknown input number'); end
else
    error('unknown inputs')
end