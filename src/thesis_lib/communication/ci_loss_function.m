function J = ci_loss_function(R1,R2,Psi,w)
% --- ci_loss_function() --------------------------------------------------
% CI loss function w.r.t. R1, R2, Psi and w.
%
% 2023-10-30 Robin Forsling

if w >= 1; J = trace(R1); return;
elseif w <= 0; error('ci_loss_function: omega out of range')
end

J = trace(inv(w*inv(R1) + (1-w)*Psi'/(Psi*R2*Psi')*Psi));