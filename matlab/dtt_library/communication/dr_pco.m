function Psi = dr_pco(R2,m)
% --- dr_pco() ------------------------------------------------------------
% Dimension-reduction (DR) based on principal component optimization (PCO).
%
% 2023-10-30 Robin Forsling

[U,D] = eig(R2);
idx_vec = get_min_idx_vec(D,m);
Psi = U(:,idx_vec)';
