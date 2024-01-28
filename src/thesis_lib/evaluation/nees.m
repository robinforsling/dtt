function [mean_nees,varargout] = nees(x,data,varargin)
% --- nees() --------------------------------------------------------------
% Computes normalized estimation error squared (NEES), the average NEES
% (ANEES), and conservativeness index (COIN).
%
% 2024-01-08 Robin Forsling

% HANDLE INPUTS
XHAT = data.xhat;
PP = data.P;

% PARAMETERS
M = length(XHAT); 
N = size(x,2);
nx = size(XHAT{1},1);
if nargin > 2; idx = varargin{1}; nx = length(idx);
else; idx = 1:nx;
end

% COMPUTE NEES AND COIN
mean_nees = NaN(1,N);
coin = NaN(1,N);
for k = 1:N
    C = zeros(nx);
    for i = 1:M
        P = PP{i}(idx,idx,k);
        Li = inv(chol(P,'lower')); 
        xtilde = XHAT{i}(idx,k)-x(idx,k);
        C = C + Li*(xtilde*xtilde')*Li';
    end
    C = C/M;
    mean_nees(k) = trace(C);
    coin(k) = max(eig(C));
end
anees = mean_nees/nx;

% COMPUTE ANEES
if nargout > 1
    varargout{1} = anees;
    varargout{2} = coin;
end