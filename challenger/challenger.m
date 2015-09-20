function [alpha, beta] = challenger(ns)
	assert(isa(ns,'double'));

	x = [53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81];
	y = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0];
	alphahat = 15.0429;
	betahat = -0.2322;
	sigmabsq = 1e-3;
	sigmaasq = 4;
	bhat = 6.0776e+06;

	alpha = zeros(ns,1);
	beta = zeros(ns,1);
	alpha(1) = alphahat+(rand*0.4-0.2);
	beta(1) = betahat+(rand*0.4-0.2);

	for i=1:ns-1
		[a,b] = ggen(alpha(i),beta(i),sigmaasq,sigmabsq);
		r = rand();
		rat = logprob(a,b,bhat,x,y) - ...
				logprob(alpha(i),beta(i),bhat,x,y);
		if log(r) < rat
			alpha(i+1) = a; 
			beta(i+1) = b;
		else
			alpha(i+1) = alpha(i); 
			beta(i+1) = beta(i);
		end
    end
end

function [alpha,beta] = ggen(a,b,sigmaasq,sigmabsq)
	alpha = a + randn()*sqrt(sigmaasq);
	beta = b + randn()*sqrt(sigmabsq);
end

function lp = logprob(a,b,bhat,x,y)
	lp = logprior(a,bhat) + logL(a,b,x,y);
end

function lp = logprior(a,bhat)
	lp = log(1/bhat) + a - exp(a)/bhat;
end

function lp = logL(a,b,x,y)
	n = length(x);
	lp = 0;
	for i=1:n
		if y(i)==1
			lp = lp + log((exp(a+b*x(i))/(1+exp(a+b*x(i)))));
		else
			lp = lp + log(1/(1+exp(a+b*x(i))));
		end
	end
end
