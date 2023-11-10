function test_message_coding
% --- test_message_coding() -----------------------------------------------
% Test the message encoding/decoding functionality.
%
% 2023-11-02 Robin Forsling


% --- PARAMETERS ---
m = 3;
nx = 6;


% --- COMPUTE PSI ---
y2 = rand(nx,1);
R1 = get_random_covariance(nx);
R2 = get_random_covariance(nx);
Q = R1^2;
S = R1+R2;
[X,E] = eig(Q,S);
idx = get_max_idx_vec(E,m);
Psi = gram_schmidt_process(X(:,idx)');
[U,D] = eig(Psi*R2*Psi');
Psi = U'*Psi;


% --- ENCODING ---
fprintf('\n\n--- PARAMETERS TO ENCODE ---:\n')
yPsi = Psi*y2
RPsi = Psi*R2*Psi'
Psi

fprintf('--- ENCODED MESSAGE ---\n')
msg = message_encoding(yPsi,RPsi,Psi)


% --- DECODING ---
fprintf('--- DECODED MESSAGE ---\n')
[yPsi_d,RPsi_d,Psi_d] = message_decoding(msg)