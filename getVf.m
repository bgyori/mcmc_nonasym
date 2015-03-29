function Vf = getVf(f)
% GETVF: Calculates the variance of the function f based on values
% f: vector or matrix of function values
	if isvector(f)
		f = f(:);
	end
	n = size(f,1);
	Vf = (1/n)*(sum(f.^2)) - ((1/n)*sum(f)).^2;
end