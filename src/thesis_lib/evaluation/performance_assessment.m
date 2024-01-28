function res = performance_assessment(x,data)
% --- performance_assessment() --------------------------------------------
% Measures for evaluation of estimator performance.
%
% 2024-01-15 Robin Forsling

% INPUTS:
XHAT = data.xhat;
PP = data.P;

% PARAMETERS:
M = length(XHAT); 
[~,N] = size(x);
nx = size(XHAT{1},1);

% POSITION COMPONENTS:
rmt.pos = zeros(1,N); rmse.pos = rmt.pos; idx = 1:2;
for k = 1:N
    r = 0; P = zeros(2);
    for i = 1:M
        xtilde = XHAT{i}(idx,k)-x(idx,k);
        r = r + xtilde'*xtilde;
        P = P + PP{i}(idx,idx,k);
    end
    rmt.pos(k) = sqrt(trace(P/M));
    rmse.pos(k) = sqrt(r/M);
end

% VELOCITY COMPONENTS:
if nx >= 4
    rmt.vel = zeros(1,N); rmse.vel = rmt.vel; idx = 3:4;
    for k = 1:N
        r = 0; P = zeros(2);
        for i = 1:M
            xtilde = XHAT{i}(idx,k)-x(idx,k);
            r = r + xtilde'*xtilde;
            P = P + PP{i}(idx,idx,k);
        end
        rmt.vel(k) = sqrt(trace(P/M));
        rmse.vel(k) = sqrt(r/M);
    end
end

% ACCELERATION COMPONENTS:
if nx >= 6
    rmt.acc = zeros(1,N); rmse.acc = rmt.acc; idx = 5:6;
    for k = 1:N
        r = 0; P = zeros(2);
        for i = 1:M
            xtilde = XHAT{i}(idx,k)-x(idx,k);
            r = r + xtilde'*xtilde;
            P = P + PP{i}(idx,idx,k);
        end
        rmt.acc(k) = sqrt(trace(P/M));
        rmse.acc(k) = sqrt(r/M);
    end
end

% SET OUTPUT RESULTS:
res.rmt = rmt;
res.rmse = rmse;

