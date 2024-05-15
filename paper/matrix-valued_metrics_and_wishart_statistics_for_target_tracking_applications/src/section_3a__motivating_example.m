function section_3a__motivating_example
% --- section_3a__motivating_example() ------------------------------------
% Motivation example of Section 3-A
%
% 2024-05-15 Robin Forsling


addpath('paper_lib/')
addpath(genpath('../../../src/thesis_lib/'))

set_latex_interpreter;


% --- SIMULATION PARAMETERS -----------------------------------------------
M = 1e5;
nk = 10;
kvec = 1:nk;


% --- MODEL PARAMETERS ----------------------------------------------------
mu = zeros(2,1);

Sigma = [8 1 ; 1 2];
P = [8 0 ; 0 2];

L = chol(P,'lower'); Li = eye(2) / L;
Xi = Li*Sigma*Li';
lambda = eig(Xi);
lambdamax = max(lambda);
lambdamin = min(lambda);


% --- PREALLOCATE ---------------------------------------------------------
nees = zeros(1,nk); nees_var = nees; anees = nees; anees_var = nees;
coin = nees; pein = nees;
Xihat = zeros(2,2,nk);
XTILDE = [];


% --- RUN MC SIMULATIONS --------------------------------------------------
for k = 1:nk 
    xtilde = sample_Gaussian_RV(mu,Sigma,M);
    stats = compute_nees_statistics(xtilde,P);
    nees(k) = stats.nees;
    nees_var(k) = stats.nees_var;
    anees(k) = stats.anees;
    anees_var(k) = stats.anees_var;
    Xihat(:,:,k) = stats.Xi;
    Lambda = eig(Xihat(:,:,k));
    coin(k) = max(Lambda);
    pein(k) = min(Lambda);

    XTILDE = [XTILDE xtilde];
end

hist_e2 = compute_histogram(XTILDE,P);
xbin = hist_e2.xbin;
pdf_samp = hist_e2.pdf_samp;
pdf0 = hist_e2.pdf0;


% --- RESULTS -------------------------------------------------------------
figure(1);clf;
clr = get_thesis_colors;

subplot(1,2,1);hold on
h1 = plot(kvec,ones(1,nk),'k--'); h1.LineWidth = 2;
hmax = plot(kvec,lambdamax*ones(1,nk),'k:');  
hmin = plot(kvec,lambdamin*ones(1,nk),'k-.');  
hanees = plot(kvec,anees,'-'); hanees.LineWidth = 2; hanees.Color = clr.darkgreen;
xlabel('$k$','interpreter','latex'); 
legend(gca,[hmax hmin hanees],'$\lambda_{\max}$','$\lambda_{\min}$','NEES/2');
box on

subplot(1,2,2);hold on
h0 = plot(xbin,pdf0,'k--'); h0.LineWidth = 2;
hsamp = plot(xbin,pdf_samp,'-'); hsamp.LineWidth = 2; hsamp.Color = clr.yellow;
xlabel('$k$','interpreter','latex');
legend(gca,[h0 hsamp],'$\chi^2$ PDF','sampled PDF')
box on

end





% --- LOCAL FUNCTIONS -----------------------------------------------------
function stats = compute_nees_statistics(xtilde,P)
    [nx,M] = size(xtilde); 
    L = chol(P,'lower'); Li = eye(nx) / L;
    Xi = Li*(xtilde*xtilde')*Li'/M;
    nu = zeros(1,M);
    for i = 1:M
        nu(i) = xtilde(:,i)'/P*xtilde(:,i);
    end
    stats.nees = mean(nu);
    stats.nees_var = var(nu);
    stats.anees = stats.nees/nx;
    stats.anees_var = var(nu/nx);
    stats.Xi = Xi;
end

function h = compute_histogram(xtilde,P)  
    [nx,N] = size(xtilde);
    nbins = 100; 

    % --- SAMPLED PDF ---
    ne2 = zeros(1,N);
    for k = 1:N
        ne2(k) = xtilde(:,k)'/P*xtilde(:,k);
    end
    emin = 0; emax = max(ne2);
    edges = linspace(emin,emax,nbins+1);
    h.xbin = edges(1:nbins) + (edges(2)-edges(1));
    h.xedge = edges;
    h.pdf_samp = histcounts(ne2,edges)/N;

    % --- TRUE PDF ---
    h.pdf0 = zeros(1,nbins);
    for i = 1:nbins
        h.pdf0(i) = chi2cdf(edges(i+1),nx) - chi2cdf(edges(i),nx);
    end
end

function X = sample_Gaussian_RV(mu,Sigma,varargin)
    if nargin > 2; M = varargin{1}; else; M = 1; end
    L = chol(Sigma,'lower');
    X = L*randn(length(mu),M) + mu;
end