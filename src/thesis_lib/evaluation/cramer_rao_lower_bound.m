function res = cramer_rao_lower_bound(X,sensor_pos,sm,pm)
% --- cramer_rao_lower_bound() --------------------------------------------------------------
% Computes Cramer-Rao lower bound (CRLB) based on input model:
%   sm = sensor model
%   pm = process model
% Both agent local CRLB and a global CRLB are computed.
%
% 2024-01-16 Robin Forsling

nsensor = length(sm);
F = pm.F;
Q = pm.Q;
nx = pm.nx;
v_max = pm.v_max;
a_max = pm.a_max;
m = sm{1}.m;
Hc = zeros(m,nx-m);
N = size(X,2);

% --- INIT ---
PGi = zeros(m);
PL = cell(nsensor,1);
for i = 1:nsensor
    xrel = X(1:2,1)-sensor_pos{i};
    H = sm{i}.J(xrel);
    % LOCAL:
    PL{i} = blkdiag(inv(H'/sm{i}.R*H),v_max^2*eye(2));
    if nx/2 >= 3; PL{i} = blkdiag(PL{i},a_max^2*eye(2)); end
    % GLOBAL:
    PGi = PGi + H'/sm{i}.R*H;
end
PG = blkdiag(inv(PGi),v_max^2*eye(2));
if nx/2 >= 3; PG = blkdiag(PG,a_max^2*eye(2)); end

% --- ALLOCATE ---
glob.P = zeros(nx,nx,N);
glob.P(:,:,1) = PG;
glob.rt.pos = zeros(1,N);
if nx >= 4; glob.rt.vel = NaN(1,N); end
if nx >= 6; glob.rt.acc = NaN(1,N); end
loc = cell(nsensor,1);
for i = 1:nsensor
    loc{i}.P = zeros(nx,nx,N);
    loc{i}.P(:,:,1) = PL{i};
    loc{i}.rt.pos = zeros(1,N);
    if nx >= 4; loc{i}.rt.vel = NaN(1,N); end
    if nx >= 6; loc{i}.rt.acc = NaN(1,N); end
end

% --- MAIN LOOP ---
for k = 2:N
    PG = F*PG*F' + Q;
    PGi = inv(PG);
    for i = 1:nsensor
        xrel = X(1:2,k)-sensor_pos{i};
        H = [sm{i}.J(xrel) Hc];
        % LOCAL:
        PL{i} = F*PL{i}*F' + Q;
        PL{i} = inv(inv(PL{i})+H'/sm{i}.R*H);
        loc{i}.P = PL{i};
        loc{i}.rt.pos(k) = sqrt(trace(PL{i}(1:2,1:2)));
        if nx >= 4; loc{i}.rt.vel(k) = sqrt(trace(PL{i}(3:4,3:4))); end
        if nx >= 6; loc{i}.rt.acc(k) = sqrt(trace(PL{i}(5:6,5:6))); end
        % GLOBAL:
        PGi = PGi + H'/sm{i}.R*H;
    end
    PG = inv(PGi);
    glob.P(:,:,k) = PG;
    glob.rt.pos(k) = sqrt(trace(PG(1:2,1:2)));
    if nx >= 4; glob.rt.vel(k) = sqrt(trace(PG(3:4,3:4))); end
    if nx >= 6; glob.rt.acc(k) = sqrt(trace(PG(5:6,5:6))); end
end

% --- SET OUTPUT RESULTS ---
res.glob = glob;
res.loc = loc;

