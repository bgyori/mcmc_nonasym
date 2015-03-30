compile;

% Parameters of Challenger model
xdata = [53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81];
ydata = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0];
alphahat = 15.0429;
alphasigma = sqrt(4);
betahat = -0.2322;
betasigma = sqrt(100e-5);
bhat= 6.0776e6;

% Estimate parameters from one run
ns = 1e7;
t0 = 0.1*ns;
[a,b] = challenger(ns+t0, xdata, ydata, alphahat, alphasigma, ...
				betahat, betasigma, bhat, 1234);
gamma = getGamma([a((t0+1:end)),b(t0+1:end)]);
eab30 = exp(a((t0+1):end)+b((t0+1):end)*30);
eab95 = exp(a((t0+1):end)+b((t0+1):end)*95);
y30 = eab30./(1+eab30);
y95 = eab95./(1+eab95);
sigma30 = getSigmaNonasym(y30',ceil(10*ns^(1/3)));
sigma95 = getSigmaNonasym(y95',ceil(10*ns^(1/3)));
V30 = getVf(y30);
V95 = getVf(y95);

challengers = load('challenger.mat');
C = 1;
ns = challengers.ns-challengers.t0;

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