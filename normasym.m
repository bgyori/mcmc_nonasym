function logp = normasym(t,n,sigma)
% NORMASYM: returns the logarithm of the error bound derived from the 
% asymptotic normal assumption
% n: number of steps (excluding burn-in)
% sigma: asymptotic or non-asymptotic variance
	nt = normcdf(t,0,sqrt(sigma/n));
	logp = zeros(size(t));
	logp(t<0) = log(nt(t<0));
	logp(t>=0) = log(1-nt(t>=0));
end