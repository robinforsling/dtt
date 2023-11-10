function out_data = example_3

options = sdpsettings('solver','mosek','verbose',0,'debug',0);

T2D = @(a) [cos(a) -sin(a) ; sin(a) cos(a)];                                % Rotation matrix
T60 = T2D(pi/3); T120 = T2D(2*pi/3);

H = [eye(2) ; eye(2) ; eye(2)];
R1 = [16 0 ; 0 1];
R2 = T60*R1*T60';
R3 = T120*R1*T120';

P = sdpvar(2);                                                              % Optimization variable
K = sdpvar(2,6);                                                            % Optimization variable
R12 = sdpvar(2,2,'full');                                                   % Uncertainty
R13 = sdpvar(2,2,'full');                                                   % Uncertainty
R23 = sdpvar(2,2,'full');                                                   % Uncertainty
R = [R1 R12 R13 ; R12' R2 R23 ; R13' R23' R3];

F = [K*H == eye(2), uncertain(R12), uncertain(R13), uncertain(R23)];        % Constraints
F = [F, R >= 0, [P K*R ; R*K' R] >= 0];                                     % Constraints (cont.)
J = trace(P);                                                               % Loss function
optimize(F, J, options);                                                    % Optimize

out_data.K = value(K);
out_data.P = value(P);