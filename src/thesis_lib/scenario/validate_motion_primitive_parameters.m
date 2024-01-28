function is_valid = validate_motion_primitive_parameters(par)
% --- validate_motion_primitive_parameters() ------------------------------
% Validates the struct par consisting of motion primitive parameters.
%
% 2023-10-30 Robin Forsling

is_valid = 1;
Ld_tol = 1e-6;

% Check fields are the same
if ~are_fieldnames_same(par,motion_primitive_parameters)
    is_valid = 0; return
end

% Check for empty fields
if isempty(par.N) || isempty(par.X) || isempty(par.head) || isempty(par.L) || isempty(par.Ld)
    is_valid = 0; return;
end

% Check for wrong N
if size(par.X,2) ~= par.N
    is_valid = 0; return;
end
if length(par.head) ~= par.N
    is_valid = 0; return;
end

% Check discretized length
L_discrete = 0;
for k = 2:par.N
    L_discrete = L_discrete + norm(par.X(1:2,k)-par.X(1:2,k-1));
end

if abs((par.Ld-L_discrete)/par.Ld) > Ld_tol
    is_valid = 0; return;
end

% Check length L
if par.Ld > par.L
    is_valid = 0; return;
end




