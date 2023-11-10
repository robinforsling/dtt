function out_data = example_5(varargin)

options = sdpsettings('solver','mosek','verbose',0,'debug',0);

rhomax = 1;                                                                 % Default rho
if nargin > 0; rhomax = varargin{1}; end
if rhomax > 1 || rhomax < 0; warning('rhomax is out of range: forcing rhomax = 0'); rhomax = 0; end                                

H1 = [1 1]/sqrt(2);                                                                                                                                                                                                                                                                                                                                                                                                              
H2 = [[1 -1]/sqrt(2) ; eye(2)];
H = [H1 ; H2];
R1 = 1;
R2 = 4*eye(3);

P = sdpvar(2);                                                              % Optimization variable                                                            
K = sdpvar(2,4);                                                            % Optimization variable                                                     
R12 = sdpvar(1,3,'full');                                                   % Uncertainty
R = [R1 R12 ; R12' R2];

F = [K*H == eye(2), uncertain(R12), [rhomax^2*R2 R12' ; R12 R1] >= 0];     % Constraints             
F = [F, R >= 0, [P K*R ; R*K' R] >= 0];                                     % Constraints (cont.)
J = trace(P);                                                               % Loss function
optimize(F, J, options);                                                    % Optimize

out_data.K = value(K);
out_data.P = value(P);
out_data.rhomax = rhomax;