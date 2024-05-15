function F = cdf_lambdamax_Wishart_gamma_approx(m,n,bvec)
% --- cdf_lambdamax_Wishart_gamma_approx() --------------------------------
% Computes an approximate CDF of the largest eigenvalues of an m x m real 
% matrix M ~ W(n,I), where W(n,I) is the real Wishart distribution of n 
% degrees of freedom and covariance parameter I. The approximation is based
% on a shifted gamma distribution.
%
% The approximation is proposed in the paper:
% M. Chiani, "Distribution of the largest eigenvalue for real Wishart and
% Gaussian random matrices and a simple approximation for the Tracy-Widom
% distribution", Journal of Multivariate Analysis, vol. 129, pp. 69–81,
% 2014.
%
% 2024-03-20 Robin Forsling


if m > n; error('m cannot be larger n...'); end


% --- TRACY-WIDOM PARAMETERS ---
mu1 = -1.206533574;                                                         % Mean of TW1. Wiki: -1.2065335745820
sigmasq1 = 1.607781034;                                                     % Variance of TW1. Wiki: 1.607781034581
s1 = 0.293464524;                                                           % Skew of TW1. Wiki: 0.29346452408

kappa = 4/s1^2;
theta = sqrt(sigmasq1)*s1/2;
rho = kappa*theta - mu1;

c1 = -1/2;                                                                  % Adjustment parameter
c2 = c1;                                                                    % Adjustment parameter
mu = (sqrt(m+c1) + sqrt(n+c2))^2;
sigma = sqrt(mu)*(1/sqrt(m+c1) + 1/sqrt(n+c2))^(1/3);


% --- APPROXIMATE CDF ---
F = zeros(size(bvec));
for i = 1:length(bvec)
    x = (bvec(i)-mu)/sigma;
    F(i) = approx_cdf(x,rho,theta,kappa);
end

end

function F = approx_cdf(x,rho,theta,kappa)
    if (x+rho) > 0; F = gammainc((x+rho)/theta,kappa);
    else; F = 0;
    end
end
