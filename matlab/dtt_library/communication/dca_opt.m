function [Ds2,varargout] = dca_opt(R2)
% --- dca_opt() -----------------------------------------------------------
% Diagonal covariance approximation using optimization based scaling. The
% optimization problem is solved using YALMIP with the MOSEK solver.
%
% 2023-10-30 Robin Forsling

nx = size(R2,1);
d2 = diag(R2);
D2 = diag(d2);
Ti = sqrtm(D2);
T = inv(Ti);
C2 = T*R2*T;                                                                % Transform to the correlation domain (to avoid numerical problems)

c = sdpvar(nx,1);                                                           % Define optimization variables
F = [diag(c)-C2 >= 0, c(:) >= 1, c(:) <= nx];                               % Constraints
options = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);       % Solver options
J = sum(c);                                                                 % J = trace
res = optimize(F,J,options);                                                % Solve the problem

s = value(c);                                                               % Scaling factors
Ds2 = Ti*diag(s)*Ti;                                                        % Transforms back to original domain
                                                              
if nargout > 1; varargout{1} = s; end
    