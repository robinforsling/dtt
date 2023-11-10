function out_var = def_association_method(varargin)
% --- def_association_method(varargin) ------------------------------------
% Definition of association method. If input is string then output is the
% corresponding association method number. If input is a discrete scalar 
% then the output is the corresponding association method string.
%
% 2023-10-30 Robin Forsling

defs = {0, 'bypass'; ...
        1, 'nn'; ...
        2, 'gnn'};

% HANDLE INPUTS
input_is_string = 0; input_is_scalar = 0;
if nargin > 0
    var = varargin{1};
    if isstring(var) || ischar(var); input_is_string = 1; 
    elseif isnumeric(var) && isscalar(var); input_is_scalar = 1; end
else
    out_var = defs; return;
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