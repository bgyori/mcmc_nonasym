compile;

nr = 1e5;
ns = 1e2;
t0 = 2000;
y30 = zeros(nr,1);
y95 = zeros(nr,1);
parfor i=1:nr
    abarray = challenger_mex(ns+t0);
	a = abarray(:,1);
    b = abarray(:,2);
    eab30 = exp(a((t0+1):end) + b((t0+1):end) * 30);
	eab95 = exp(a((t0+1):end) + b((t0+1):end) * 95);
    y30(i) = mean(eab30 ./ (1+eab30));
	y95(i) = mean(eab95 ./ (1+eab95));
end
save('challenger_out.mat');
