function [x,varargout] = inv_cdf_lambdamin_Wishart(m,n,P,varargin)
% --- inv_cdf_lambdamin_Wishart() -----------------------------------------
% Computes the inverse of the CDF of the smallest eigenvalue of a matrix 
% M ~ W(n,I), where W(n,I) is the real Wishart distribution of n degrees of
% freedom and with covariance parameter I.
%
% 2024-03-21 Robin Forsling


% --- PARAMETERS ---
m_threshold = 9;
n_threshold = 60;
P_threshold = 1e-5;
max_iter = 100;

if P >= 1; error('P should lie in the semiopen interval [0,1)'); end
if m > n; error('m cannot be larger n...'); end
if nargin > 3 && ~isempty(varargin{1}); P_threshold = varargin{1}; end
if nargin > 4 && ~isempty(varargin{2}); max_iter = varargin{2}; end
if nargin > 5 && ~isempty(varargin{3}); n_threshold = varargin{3}; end


% --- SELECT CDF ---
if m > m_threshold || n > n_threshold; CDF = @cdf_lambdamin_Wishart_gamma_approx;
else; CDF = @cdf_lambdamin_Wishart; 
end


% --- INITIAL GUESS ---
xmin = 0; xmax = n;
while P > CDF(m,n,xmax)
    xmin = xmax; xmax = 2*xmax;
end


% --- REFINED SEARCH ---
cnt_iter = 0;
F = CDF(m,n,xmax);
while abs(F-P) > P_threshold && cnt_iter < max_iter
    xmid = 0.5*(xmin+xmax);
    F = CDF(m,n,xmid);
    if P > F; xmin = xmid;
    else; xmax = xmid;
    end
    cnt_iter = cnt_iter + 1;
end

x = xmid;
if nargout > 1; varargout{1} = cdf_lambdamin_Wishart(m,n,x); end

