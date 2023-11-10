function Ds2 = dca_eig(R2)
% --- dca_eig() -----------------------------------------------------------
% Diagonal covariance approximation using eigenvalue based scaling.
%
% 2023-10-30 Robin Forsling

D2 = diag(diag(R2));
T = inv(sqrt(D2));
s = max(eig(T*R2*T));
Ds2 = s*D2;