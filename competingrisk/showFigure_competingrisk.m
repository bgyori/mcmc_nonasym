addpath ..
compile;

% Estimate parameters from one run
ns = 1e5;
t0 = 20000;
[out,outend,data] = competingrisk_mex(ns+t0);
gamma = getGamma(out(t0+1:end,:));

mxdiv = 200;
mxt = 20;
dt = mxt/mxdiv;
t = linspace(mxt/mxdiv,mxt,mxdiv);
a = 1;

[r, nstages] = size(out);
f = zeros(ns,1);
for i = 1:ns
	b1 = out(t0+i, nstages-3);
	sig1 = out(t0+i, nstages-2);
	b2 = out(t0+i, nstages-1);
	sig2 = out(t0+i, nstages);      


	h1 = b1*t.^(b1-1)./(sig1^b1*(1+(t/sig1).^b1));
	h2 = b2*t.^(b2-1)./(sig2^b2*(1+(t/sig2).^b2));
	s1 = 1./(1+(t/sig1).^b1);
	s2 = 1./(1+(t/sig2).^b2);

	pd1 = (h1+h2).*s1.*s2;
	pd0 = s1.*s2;
	spd1 = sum(pd1);
	spd0 = sum(pd0);
	psum = spd1+spd0;

	p = (pd1+pd0)/psum;
	Et = sum(p.*t);
	f(i) = (Et>a);
end

sigma = getSigmaNonasym(f, ceil(10*ns^(1/3)));
V = getVf(f);

competing = load('competingrisk.mat');
C = 1;
ns = competing.ns;

[t,logp] = logTail(competing.probtla,100);
logp_bern = bernsteintail(t,ns,sigma,V,gamma,C);
logp_cheb = chebyshevtail(t,ns,sigma,V,gamma);
logp_norm = normasym(t,ns,sigma);

fprintf('sigma^2 = %.3e\tVf = %.3e\tgamma = %.3e\n',sigma,V,gamma);
logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm);