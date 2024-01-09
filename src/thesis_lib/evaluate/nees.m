function [mean_nees,varargout] = nees(x,data,varargin)
% --- nees() --------------------------------------------------------------
% Computes normalized estimation error squared (NEES), the average NEES
% (ANEES), and conservativeness index (COIN).
%
% 2024-01-08 Robin Forsling

% HANDLE INPUTS
if isstruct(data)
    if isfield(data,'xhat'); XHAT = data.xhat; 
    elseif isfield(data,'XHAT'); XHAT = data.XHAT; 
    elseif isfield(data,'XX'); XHAT = data.XX; 
    else; error('invalid input, cannot retrieve xhat')
    end
    if nargin > 1
        if isfield(data,'P'); PP = data.P; 
        elseif isfield(data,'PP'); PP = data.PP; 
        else; error('invalid input, cannot retrieve P')
        end
    end
else
    error('invalid input')
end
if ~(iscell(XHAT) && iscell(PP)); error('invalid input, struct fields must be cells'); end
if length(XHAT) ~= length(PP); error('invalid input, lengths of cell arrays must be equal'); end

% PARAMETERS
M = length(XHAT); 
N = size(x,2);
nx = size(x,1);
if nargin > 2; idx = varargin{1}; nx = length(idx);
else; idx = 1:nx;
end
x = x(idx,1:N);

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



% % COMPUTE MEAN NEES
% mean_nees = NaN(1,N);
% for k = 1:N
%     s = 0;
%     for i = 1:M
%         dx = XHAT{i}(idx,k)-x(idx,k);
%         P = PP{i}(idx,idx,k);
%         s = s + dx'/P*dx;
%     end
%     mean_nees(k) = s/(M);
% end
% 
% % COMPUTE ANEES
% if nargout > 1
%     varargout{1} = mean_nees/nx;
% end