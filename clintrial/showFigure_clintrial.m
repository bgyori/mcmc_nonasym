addpath ..
compile;

% Single run to calculate parameters
ns = 1e5;
t0 = 600;
f = clintrial_mex(ns+t0);
[~, allvariables] = clintrial_mex(ns+t0);
gamma = getGamma(allvariables(t0+1:end-1,:));
sigma = getSigmaNonasym(f(t0+1:end),ceil(10*ns^(1/3)));
V = getVf(f(t0+1:end));

% Load results
clin = load('clintrial_out.mat');
C = 1;
ns = clin.ns;

% Plot error bounds
[t,logp] = logTail(clin.avarray,100);
logp_bern = bernstein_nonrev(t,ns,sigma,V,gamma,C);
logp_cheb = chebyshev_nonrev(t,ns,sigma,V,gamma);
logp_norm = normasym(t,ns,sigma);

fprintf('sigma^2 = %.3e\tVf = %.3e\tgammaps = %.3e\n',sigma,V,gamma);
logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm);