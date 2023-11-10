function J = kf_loss_function(R1,R2,Psi)
% --- kf_loss_function() --------------------------------------------------
% KF loss function w.r.t. R1, R2 and Psi.
%
% 2023-10-30 Robin Forsling

J = trace(R1 - R1*Psi'/(Psi*(R1+R2)*Psi')*Psi*R1);