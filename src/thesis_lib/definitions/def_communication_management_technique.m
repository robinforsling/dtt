function out_var = def_communication_management_technique(varargin)
% --- def_communication_management_technique(varargin) --------------------
% Definition of communication management technique, e.g., GEVO and DCA
%
% 2023-10-30 Robin Forsling

defs = {0, 'no-communication'; ...                                          % No estimates communicated
        1, 'full'; ...                                                      % Full estimates communicated    
        2, 'pco'; ...                                                       % Principle component optimization (PCO)
        3, 'paro'; ...                                                      % Principle axis restricted optimization (PARO, cf. FUSION2020)
        4, 'gevo'; ...                                                      % Generalized eigenvalue optimization (GEVO, cf. FUSION2022)
        5, 'dca-eig'; ...                                                   % Diagonal covariance approximation (DCA) - eigenvalue based scaling
        6, 'dca-opt'; ...                                                   % Diagonal covariance approximation - optimization based scaling
        7, 'dca-dom'; ...                                                   % Diagonal covariance approximation - diagonal-dominance based scaling
        8, 'dca-hyp'};                                                      % Diagonal covariance approximation - hyperrectangle enclosing

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