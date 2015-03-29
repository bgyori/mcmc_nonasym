function logp = bernstein_nonrev(t,n,sigma,Vf,gamma,C)
% BERNSTEIN: returns the logarithm of the Bernstein-type error bound for
% non-reversible Markov chains. 
% n: number of steps (excluding burn-in)
% sigma: non-asymptotic variance
% Vf: variance of the function
% gamma: spectral gap
% C: bound on the value of the function
	logp = -(n-1/gamma)*gamma*t.^2 ./ (8*Vf + 20*C*abs(t));
end