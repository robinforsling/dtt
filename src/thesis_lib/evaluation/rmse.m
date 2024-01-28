function [r,varargout] = rmse(x,data)
% --- rmse() --------------------------------------------------------------
% Computes root mean squared error (RMSE) and root mean trace (RMT).
%
% 2023-10-30 Robin Forsling

% HANDLE INPUTS
XHAT = data.xhat;
PP = data.P;

% PARAMETERS
M = length(XHAT); 
N = size(x,2);
nx = size(x,1);
p_idx = 1:2;
v_idx = 3:4;
a_idx = 5:6;

% COMPUTE RMSE
r = NaN(1,N);
for k = 1:N
    s = 0;
    for i = 1:M
        dx = XHAT{i}(p_idx,k)-x(p_idx,k);
        s = s + dx'*dx;
    end
    r(k) = sqrt(s/M);
end

% COMPUTE SQRT MEAN TRACE P
if nargout > 1
    stp = NaN(1,N);
    for k = 1:N
        s = 0;
        for i = 1:M
            P = PP{i}(p_idx,p_idx,k);
            s = s + trace(P);
        end
        stp(k) = sqrt(s/M);
    end
    varargout{1} = stp;
end


