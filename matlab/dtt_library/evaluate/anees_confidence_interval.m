function conf_int = anees_confidence_interval(n,M)
% --- anees_confidence_interval() -----------------------------------------
% Constructs confidence intervals for ANEES.
%
% 2023-10-30 Robin Forsling

d = 9*n*M;
b = 1.96*sqrt(2/d);
conf_int.lvl95 = [1 - 2/d - b , 1 - 2/d + b].^3;

b = 2.576*sqrt(2/d);
conf_int.lvl99 = [1 - 2/d - b , 1 - 2/d + b].^3;

b = 3.2905*sqrt(2/d);
conf_int.lvl999 = [1 - 2/d - b , 1 - 2/d + b].^3;