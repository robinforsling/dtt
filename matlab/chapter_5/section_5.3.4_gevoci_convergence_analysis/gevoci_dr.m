function k = gevoci_dr(R1,R2,m,J_thres)
% --- gevoci_dr() ---------------------------------------------------------
% GEVO for covariance intersection
%
% 2023-10-30 Robin Forsling

max_iter = 50;
wmin = 0.001;
wmax = 1;
I1 = inv(R1);

% 1. Initialize
w = 0.5;
J = Inf;

for k = 1:max_iter

    J_prev = J;

    % 2. Solve for X
    A = R1/w; B = R2/(1-w);
    Q = A^2; S = A+B;
    [X,D] = eig(Q,S);
    idx = get_max_idx_vec(D,m);
    Psi = X(:,idx)';
    RPsi = Psi*R2*Psi';
    
    % 3. Solve for w
    f = @(w) trace( inv(w*I1 + (1-w)*Psi'/RPsi*Psi) );
    w = fminbnd(f,wmin,wmax);

    J = loss_function(I1,RPsi,Psi,w);

    if abs(J_prev-J)/J < J_thres; break; end
end

end

function J = loss_function(I1,RPsi,Psi,w)
    J = trace( inv(w*I1+(1-w)*Psi'/RPsi*Psi) );
end