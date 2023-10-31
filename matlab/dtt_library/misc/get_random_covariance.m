function R = get_random_covariance(varargin)
% --- get_random_covariance() ---------------------------------------------
% Samples random covariance.
%
% 2023-10-30 Robin Forsling

if nargin > 0; n = varargin{1}; else; n = 2; end

r = randn(n);
R = r*r';