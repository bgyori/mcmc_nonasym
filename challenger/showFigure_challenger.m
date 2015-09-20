addpath ..
compile;

% Estimate parameters from one run
ns = 1e5;
t0 = 2000;
[a, b] = challenger_mex(ns+t0);
gamma = getGamma([a(t0+1:end,:), b(t0+1:end,:)]);
eab30 = exp(a((t0+1):end)+b((t0+1):end)*30);
eab95 = exp(a((t0+1):end)+b((t0+1):end)*95);
y30 = eab30./(1+eab30);
y95 = eab95./(1+eab95);
sigma30 = getSigmaNonasym(y30',ceil(10*ns^(1/3)));
sigma95 = getSigmaNonasym(y95',ceil(10*ns^(1/3)));
V30 = getVf(y30);
V95 = getVf(y95);

challengers = load('challenger_out.mat');
C = 1;
ns = challengers.ns - challengers.t0;

[t,logp] = logTail(challengers.y30,100);
logp_bern = bernsteintail(t,ns,sigma30,V30,gamma,C);
logp_cheb = chebyshevtail(t,ns,sigma30,V30,gamma);
logp_norm = normasym(t,ns,sigma30);

fprintf('sigma^2 = %.3e\tVf = %.3e\tgamma = %.3e\n',sigma30,V30,gamma);
logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm);


[t,logp] = logTail(challengers.y95,100);
logp_bern = bernsteintail(t,ns,sigma95,V95,gamma,C);
logp_cheb = chebyshevtail(t,ns,sigma95,V95,gamma);
logp_norm = normasym(t,ns,sigma95);

fprintf('sigma^2 = %.3e\tVf = %.3e\tgamma = %.3e\n',sigma95,V95,gamma);
logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm);