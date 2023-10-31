function [xhat,P] = subtract_estimate(y1,R1,y2,R2)
% --- subtract_estimate() -------------------------------------------------
% Subtracts estimate (y2,R2) from estimate (y1,R2).
%
% 2023-10-30 Robin Forsling

I1 = inv(R1); I2 = inv(R2);
P = inv(I1-I2);
xhat = P*(I1*y1-I2*y2);

if min(eig(P)) < 0; warning('subtract_estimate: P is not PSD'); disp(min(eig(P))); end
