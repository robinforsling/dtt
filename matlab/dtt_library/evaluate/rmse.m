function [r,varargout] = rmse(x,data,varargin)
% --- rmse() --------------------------------------------------------------
% Computes root mean squared error (RMSE).
%
% 2023-10-30 Robin Forsling


% HANDLE INPUTS
if isstruct(data)
    if isfield(data,'xhat'); XHAT = data.xhat; 
    elseif isfield(data,'XHAT'); XHAT = data.XHAT; 
    elseif isfield(data,'XX'); XHAT = data.XX; 
    else; error('invalid input, cannot retrieve xhat...')
    end
    if nargin > 1
        if isfield(data,'P'); PP = data.P; 
        elseif isfield(data,'PP'); PP = data.PP; 
        else; error('invalid input, cannot retrieve P...')
        end
    end
else
    error(func_name,': invalid input')
end
if ~(iscell(XHAT) && iscell(PP)); error('invalid input, struct fields must be cells'); end
if length(XHAT) ~= length(PP); error('invalid input, lengths of cell arrays must be equal'); end

% PARAMETERS
M = length(XHAT); 
N = size(x,2);
d = size(x,1);
r_idx = 1:d;
if nargin > 2; p_idx = varargin{1};
else; p_idx = 1:d;
end
x = x(r_idx,1:N);

% COMPUTE RMSE
r = NaN(1,N);
for k = 1:N
    s = 0;
    for i = 1:M
        dx = XHAT{i}(r_idx,k)-x(r_idx,k);
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


