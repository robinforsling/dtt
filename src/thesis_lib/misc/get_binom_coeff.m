function N = get_binom_coeff(n,k)
%
%
%
% 2023-10-30 Robin Forsling

N = factorial(n)/(factorial(k)*factorial(n-k));