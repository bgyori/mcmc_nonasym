function gamma = getGamma(f)
% GETGAMMA: calculates the spectral gap of the Markov chain based
% on one or more outputs in the matrix f
% f: n x nf matrix where n is the number of steps and nf is the number of 
% output functions calculated
	if isvector(f)
		f = f(:);
	end
	
	[n,nf] = size(f);
	
	Vf = getVf(f);
	
	i = 1;
	eta(1) = 1;
	while true
		for j=1:nf
			g(i,j) = getGammaC(f(:,j),eta(i),eta(i));
			if Vf(j)>0
                gammas(i,j) = 1-(g(i,j)/Vf(j)).^(1/eta(i));
            else
                gammas(i,j) = 1;
            end
		end
		ming(i) = min(gammas(i,:));
		if i>1
			if ming(i)>=ming(i-1)
				break;
			end
		end
		eta(i+1) = ceil(log(n*ming(i))/(4*log(1/(1-ming(i)))));
		i = i + 1;
	end
	gamma = ming(end);
end
