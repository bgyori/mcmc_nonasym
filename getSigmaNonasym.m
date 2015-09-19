function sigma = getSigmaNonasym(f,k)
% GETSIGMANONASYM: calculate the non-asymptotic variance based on function
% output
% f: vector containing function output from n steps of a Markov chain
% k: constant usually set to 10*n^(1/3)
	if ~isvector(f)
		error('f must be a vector')
	end
	f = f(:);
	n = length(f);
    f = f - mean(f);

    kp1sum = zeros(n-k,1);
    kp1sum(1) = sum(f(1:k+1));

    for i = 2:n-k
        kp1sum(i) = kp1sum(i-1) - f(i-1) + f(i+k);
    end

    sigma1 = sum((2*kp1sum - f(1:n-k)) .* f(1:n-k))/(n-k);
    firstav = sum(f(1:(n-k)))/(n-k);
    lastav = firstav - (sum(f(1:k) - f(n-k+1:n)))/(n-k);
    sigma = sigma1 - (2*k+1)/2*(firstav^2 + lastav^2);
    sigma = sigma*(n-k)/(n-3*k-1);
end
