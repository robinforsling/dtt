function res = credibility_assessment(x,data)
% --- credibility_assessment() --------------------------------------------
% Measures for evaluation of estimator credibility and conservativeness.
%
% 2024-01-15 Robin Forsling

% INPUTS:
XHAT = data.xhat;
PP = data.P;

% PARAMETERS:
M = length(XHAT); 
N = size(x,2);
nx = size(XHAT{1},1);
idx = 1:nx;

% COMPUTE STATISTICS:
Phat = cell(1,N); Sigmahat = Phat; Gammahat = Phat;                                            
pein = zeros(1,N); coin = pein; nees = pein; 
Y = zeros(M,N);                                                             % Normalized error vector norm
for k = 1:N
    G = zeros(nx); R = G; S = G;
    for i = 1:M
        P = PP{i}(:,:,k);
        Li = inv(chol(P,'lower')); 
        xtilde = XHAT{i}(:,k)-x(idx,k);
        ytilde = Li*xtilde;
        R = R + P;
        S = S + xtilde*xtilde';
        G = G + ytilde*ytilde';
        Y(i,k) = norm(ytilde);
    end
    Phat{k} = R/M;
    Sigmahat{k} = S/M;
    Gammahat{k} = G/M;
    lambda = eig(Gammahat{k});
    pein(k) = min(lambda);
    coin(k) = max(lambda);
    nees(k) = sum(lambda);
end
norm_of_ne = compute_histogram(Y);
norm_of_ne2 = compute_histogram(Y.^2);

% SET OUTPUT RESULTS:
res.Phat = Phat;
res.Sigmahat = Sigmahat;
res.Gammahat = Gammahat;
res.credibility_interval = [pein ; coin];
res.pein = pein;
res.coin = coin;
res.nees = nees;
res.anees = nees/nx;
res.nne = norm_of_ne;
res.nne2 = norm_of_ne2;

end


function h = compute_histogram(x)
    [M,N] = size(x);
    nbins = min([25,max([20,floor(M/20)])]);
    emin = 0; emax = max(x,[],'all');
    edges = linspace(emin,emax,nbins+1);
    h.xbin = edges(1:nbins) + (edges(2)-edges(1));
    h.xedge = edges;
    h.pdf = zeros(nbins,N);
    for k = 1:N
        h.pdf(:,k) = histcounts(x(:,k),edges)/M;
    end
    h.pdf_mean = sum(h.pdf,2)/N;
end

