function J = bsc_loss_function(R1,R2,R12,Psi)
% --- bsc_loss_function() -------------------------------------------------
% Bar-Shalom-Campo loss function w.r.t. R1, R2, R12 and Psi.
%
% 2023-10-30 Robin Forsling

A = Psi*(R1-R12');
B = Psi*(R1+R2-R12-R12')*Psi';
J = trace(R1 - A'/B*A);