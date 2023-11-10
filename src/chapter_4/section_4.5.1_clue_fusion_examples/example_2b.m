function out_data = example_2b

options = sdpsettings('solver','mosek','verbose',0,'debug',0);

H = [eye(2) ; eye(2)];
R1 = [1 0 ; 0 4];
R2 = [4 0 ; 0 1];

P = sdpvar(2);                                                              % Optimization variable
K = sdpvar(2,4);                                                            % Optimization variable
Gi = sdpvar(2);                                                             % Uncertainty, Gamma^(-1)

R = [R1 R1*Gi*R2 ; R2*Gi*R1 R2];

F = [K*H == eye(2), R >= 0, uncertain(Gi), inv(R1) >= Gi, inv(R2) >= Gi];   % Constraints
F = [F, [P K*R ; R*K' R] >= 0];                                             % Constraints (cont.)
J = trace(P);                                                               % Loss function
optimize(F, J, options);                                                    % Optimize

out_data.K = value(K);
out_data.P = value(P);