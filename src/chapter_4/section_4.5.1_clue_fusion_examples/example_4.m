function out_data = example_4

options = sdpsettings('solver','mosek','verbose',0,'debug',0);

H = [eye(2) ; eye(2)];
R1 = [5 1 ; 1 1];
R2 = [1 -1 ; -1 5];
Q = [R1 [1 0.5 ; 0.5 1] ; [1 0.5 ; 0.5 1]' R2];
S = [R1 [-1 0.5 ; 0.5 -1] ; [-1 0.5 ; 0.5 -1]' R2];

P = sdpvar(2);                                                              % Optimization variable
K = sdpvar(2,4);                                                            % Optimization variable

F = [K*H == eye(2), [P K*Q ; Q*K' Q] >= 0, [P K*S ; S*K' S] >= 0];          % Constraints
J = trace(P);                                                               % Loss function
optimize(F, J, options);                                                    % Optimize

out_data.K = value(K);
out_data.P = value(P);