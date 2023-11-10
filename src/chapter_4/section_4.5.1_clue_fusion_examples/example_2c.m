function out_data = example_2c

options = sdpsettings('solver','mosek','verbose',0,'debug',0);

H = [eye(2) ; eye(2)];
R1 = [1 0 ; 0 4];
R2 = [4 0 ; 0 1];

P = sdpvar(2);                                                              % Optimization variable
K = sdpvar(2,4);                                                            % Optimization variable
R12 = diag(sdpvar(2,1,'full'));                                             % Uncertainty
R = [R1 R12 ; R12' R2];

F = [K*H == eye(2), R >= 0, uncertain(R12), [P K*R ; R*K' R] >= 0];         % Constraints
J = trace(P);                                                               % Loss function
optimize(F, J, options);                                                    % Optimize

out_data.K = value(K);
out_data.P = value(P);