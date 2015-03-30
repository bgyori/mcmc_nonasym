addpath ..
compile

for i=1:4
	rndwalks = load(sprintf('output/rndwalk_mh%d.mat',i));
	d = length(rndwalks.sigmap);
	% Estimate parameters from one run
	n = 1e7;
	seed = 123;

	[f,xa] = rndwalk_mh(n,rndwalks.sigmap,1,2.5,seed);
	gamma = getGamma(xa);
	sigma = getSigmaNonasym(f,ceil(10*n^(1/3)));
	Vf = getVf(f);
	
	C = 1;
	n = rndwalks.ns;
	% Simulation
	[t,logp] = logTail(rndwalks.m);
	% Bernstein inequality
	logp_bern = bernsteintail(t,n,sigma,Vf,gamma,C);
	% Chebyshev inequality
	logp_cheb = chebyshevtail(t,n,sigma,Vf,gamma);
	% Normal approximation
	logp_norm = normasym(t,n,sigma);
	
	fprintf('sigma^2 = %.3e\tVf = %.3e\tgamma = %.3e\n',sigma,Vf,gamma)
	
	% Plot
	if i==1
		show_legend = true;
	else
		show_legend = false;
	end
	logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm,show_legend);
end