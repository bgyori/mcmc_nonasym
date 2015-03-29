function logp = chebyshev_nonrev(t,n,sigma,Vf,gamma)
% CHEBYSHEV_NONREV: returns the logarithm of the Chebishev-type error bound for
% non-reversible Markov chains. 
% n: number of steps (excluding burn-in)
% sigma: non-asymptotic variance
% Vf: variance of the function
% gamma: spectral gap
	p = (sigma + 16*Vf/(n*gamma^2))./(n*(t.^2));
	logp = log(min(p,1));
end