for i=1:4
	rndwalks = load(sprintf('output/cutoff3.5/rndwalk_mh%d.mat',i));
	d = length(rndwalks.sigmap);
	% Estimate parameters from one run
	n = 1e7;
	[f,xa] = rndwalk_mh(n,1,d,rndwalks.sigmap,1,3.5,ceil(1000*rand));
	gamma = getGamma(xa);
	sigma = getSigmaNonasym(f,ceil(10*n^(1/3)));
	Vf = getVf(f);
	%----
	
	C = 1;
	n = rndwalks.ns;
	% Simulation
	[t,logp] = logTail(rndwalks.m,200);
	% Bernstein inequality
	logp_bern = bernstein(t,n,sigma,Vf,gamma,C);
	% Chebyshev inequality
	logp_cheb = chebyshev(t,n,sigma,Vf,gamma);
	% Normal approximation
	logp_normal = normasym(t,n,sigma);
	
	fprintf('sigma^2 = %.3e\tVf = %.3e\tgamma = %.3e\n',sigma,Vf,gamma)
	
	% Plot
	figure; hold on;
	title(sprintf('Trial %d',i));
	set(0,'DefaultAxesColorOrder',[0 0 0]);
	set(0,'DefaultLineLineWidth',1.5);
	plot(t,logp,'-');
	plot(t,logp_cheb,'-');
	plot(t(1:5:end),logp_cheb(1:5:end),'x','MarkerSize',5);
	plot(t,logp_bern,'--');
	plot(t,logp_normal,'-.');
	ylim([min(logp(~isinf(logp)))*1.5,0]);
	xlabel('$$t$$','interpreter','latex','FontSize',15);
	ylabel('$$\widehat{L}(t)$$','interpreter','latex','FontSize',15);
	hold off;
end