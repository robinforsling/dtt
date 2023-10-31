function d2 = mahalanobis_distance(y1,R1,y2,R2)
% --- mahalanobis_distance() ----------------------------------------------
% Computes the Mahalanobis distance between random vectors y1 and y2.
%
% 2023-10-30 Robin Forsling

ytilde = y1-y2;
S = R1+R2;
d2 = ytilde'/S*ytilde;