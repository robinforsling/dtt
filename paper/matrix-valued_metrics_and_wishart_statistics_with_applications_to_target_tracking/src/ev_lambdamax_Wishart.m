function mu = ev_lambdamax_Wishart(m,n,varargin)
% --- ev_lambdamax_Wishart() ----------------------------------------------
% Computes the expected value of the largest eigenvalue of M ~ W(n,I), 
% where W(n,I) is the Wishart distribution with n degrees of freedom and
% covariance parameres I. The expected value is based on numerical
% integration.
%
% 2024-03-21 Robin Forsling


if m > n; error('m cannot be larger n...'); end


% --- PARAMETERS ---
m_threshold = 9;
n_threshold = 50;
d = 1e-3;
t = 1e-7;

if nargin > 2 && ~isempty(varargin{1}); d = varargin{1}; end
if nargin > 3 && ~isempty(varargin{2}); t = varargin{2}; end
if nargin > 4 && ~isempty(varargin{3}); n_threshold = varargin{3}; end


% --- SELECT CDF ---
if m > m_threshold || n > n_threshold; CDF = @cdf_lambdamax_Wishart_gamma_approx;
else; CDF = @cdf_lambdamax_Wishart; 
end


% --- EXPECTED VALUE ---
F = 0;
x = 0;
mu = 0;
while (1-F) > t
    F = CDF(m,n,x);
    mu = mu + d*(1-F);
    x = x + d;
end
