function ci = conf_int_lambda_Wishart(m,n,P)
% --- conf_int_lambda_Wishart() -------------------------------------------
% Returns a 100*P % confidence interval ci for the eigenvalues of a matrix 
% M ~ W(n,I), where W(n,I) is the real Wishart distribution of n degrees of
% freedom and with covariance parameter I.
%
% 2024-03-21 Robin Forsling


if m > n; error('m cannot be larger n...'); end


% --- COMPUTE CONFIDENCE INTERVAL ---
alpha = 1-P;
a = inv_cdf_lambdamin_Wishart(m,n,alpha/2);
b = inv_cdf_lambdamax_Wishart(m,n,1-alpha/2);
ci = [a b];
