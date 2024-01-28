function [PCRLB,varargout] = crlb(X,sensor_pos,sm,pm)
% --- crlb() --------------------------------------------------------------
% Computes Cramer-Rao lower bound based on input model:
%   sm = sensor model
%   pm = process model
%
% 2023-10-30 Robin Forsling

nsensors = length(sm);
F = pm.F;
Q = pm.Q;
nx = pm.nx;
d = pm.ncoord;
v_max = pm.v_max;
a_max = pm.a_max;
m = sm{1}.m;
Hc = zeros(m,nx-m);

N = size(X,2);


% --- INIT ---
Pi = zeros(m);
for i = 1:nsensors
    xrel = X(1:d,1)-sensor_pos{i};
    H = sm{i}.J(xrel);
    Pi = Pi + H'/sm{i}.R*H;
end
P = blkdiag(inv(Pi),v_max^2*eye(d));
if nx/d == 3; P = blkdiag(P,a_max^2*eye(d)); end


% --- MAIN LOOP ---
PCRLB = zeros(nx,nx,N);
PCRLB(:,:,1) = P;
for k = 2:N
    P = F*P*F' + Q;
    Pi = inv(P);
    for i = 1:nsensors
        xrel = X(1:d,k)-sensor_pos{i};
        H = [sm{i}.J(xrel) Hc];
        Pi = Pi + H'/sm{i}.R*H;
    end
    P = inv(Pi);
    PCRLB(:,:,k) = P;
end

if nargout > 1
    pcrlb = zeros(1,N);
    for k = 1:N
        pcrlb(k) = sqrt(trace(PCRLB(1:d,1:d,k)));
    end
    varargout{1} = pcrlb;
end
