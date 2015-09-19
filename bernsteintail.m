function logp = bernstein(t,n,sigma,Vf,gamma,C)
% BERNSTEIN: returns the logarithm of the Bernstein-type error bound for
% reversible Markov chains. 
% n: number of steps (excluding burn-in)
% sigma: non-asymptotic variance
% Vf: variance of the function
% gamma: spectral gap
% C: bound on the value of the function
	logp = -n*t.^2 ./ (2*(sigma+0.8*Vf)+10*C*abs(t)/gamma);
end