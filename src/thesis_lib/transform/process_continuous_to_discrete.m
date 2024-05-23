function [Fd,Qd] = process_continuous_to_discrete(Fc,Qc,dt)
% --- process_continuous_to_discrete() ------------------------------------
% Convert continuous time process to discrete time. Uses Van Loan's method 
% to work out the transformation. This produces an exact solution.
%
% Thanks to Simon Julier for sharing this code in May 2024.
%
% 2024-05-23 Robin Forsling


n = size(Fc,1);

bigA = zeros(2*n);

bigA(1:n,1:n) = -Fc*dt;
bigA(1:n, n+1:end) = Qc*dt;
bigA(n+1:end, n+1:end) = Fc'*dt;

bigB = expm(bigA);

Fd = bigB(n+1:end,n+1:end)';
Qd = Fd*bigB(1:n,n+1:end);