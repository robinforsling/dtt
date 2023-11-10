function out_var = def_fusion_method(varargin)
% --- def_fusion_method(varargin) -----------------------------------------
% Definition of fusion methods. If input is string then output is the
% corresponding fusion method number. If input is a discrete scalar then 
% the output is the corresponding fusion method string.
%
% 2023-10-30 Robin Forsling

defs = {0, 'lkf'; ...
        1, 'kf'; ...
        2, 'gimf'; ...
        3, 'bsc'; ...
        4, 'ci'; ...
        5, 'ici'; ...
        6, 'le'; ...
        7, 'ci-dca-hyp'};

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