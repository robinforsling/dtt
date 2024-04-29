function F = cdf_lambdamin_Wishart(m,n,avec)
% --- cdf_lambdamin_Wishart() ---------------------------------------------
% Computes the CDF at b of the smallest eigenvalues of an m x m real matrix 
% M ~ W(n,I), where W(n,I) is the real Wishart distribution of n degrees of 
% freedom and covariance parameter I, and a in avec. Uses log to avoid
% numerical issues.
%
% The algorithm corresponds to Algorithm 1, with b = Inf, of the paper:
% M. Chiani, "On the probability that all eigenvalues of Gaussian, Wishart, 
% and double Wishart random matrices lie within an interval", IEEE 
% Transactions on Information Theory, vol. 63, no. 76, pp. 4521-4531, 2017.
%
% 2024-03-21 Robin Forsling


if m > n; error('m cannot be larger n...'); end


% --- DEFINITIONS ---
r = @(z,a) 1-gammainc(a,z); 
glog = @(z,t) z*log(t)-t;
Gammalog = @(z) gammaln(z);
alpha = (n-m-1)/2;
Klog = normalization_constant(alpha,m,n);


% --- COMPUTE CDF ---
F = zeros(size(avec));
for k = 1:length(avec)
    a = avec(k);

    A = zeros(m);
    for i = 1:m-1
        alphai = alpha + i;
        for j = i:m-1
            alphaj = alpha + j;
            alphaij = alphai + alphaj;

            L1 = (1-alphaij)*log(2) + Gammalog(alphaij) - Gammalog(alphaj+1) - Gammalog(alphai);
            L2 = glog(alphaj,a/2) - Gammalog(alphaj+1);
            A(i,j+1) = A(i,j) + r(alphaij,a)*exp(L1) - r(alphai,a/2)*exp(L2);
        end
    end
    
    if is_odd(m)
        acol = zeros(m,1);
        for i = 1:m
            acol(i) = r(alpha+i,a/2);
        end
        A = [A acol ; zeros(1,m+1)];
    end
    
    A = A - A';
    logdetA = sum(log(eig(A)));
    F(k) = 1 - exp(Klog + (1/2)*logdetA); 
end

end

function Klog = normalization_constant(alpha,m,n)
    Klog = (m^2/2) * log(pi) - Gammamlog(m,n/2) - Gammamlog(m,m/2);
    Klog = Klog + sum(gammaln(alpha+(1:m)));
end

function val = Gammamlog(m,z)
    val = (m*(m-1)/4) * log(pi) + sum(gammaln(z-((1:m)-1)/2));
end

function val = is_odd(m)
    if mod(m,2) == 0; val = false; else; val = true; end
end