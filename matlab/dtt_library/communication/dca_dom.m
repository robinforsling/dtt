function Ds2 = dca_dom(R2)
% --- dca_dom() -----------------------------------------------------------
% Diagonal covariance approximation using diagonal-dominance based scaling. 
%
% 2023-10-30 Robin Forsling

nx = size(R2,1);
Ds2 = zeros(nx);

for i = 1:nx; Ds2(i,i) = sum(abs(R2(i,:))); end