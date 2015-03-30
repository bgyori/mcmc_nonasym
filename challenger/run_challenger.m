compile;

xdata = [53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81];
ydata = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0];
alphahat = 15.0429;
alphasigma = sqrt(4);
betahat = -0.2322;
betasigma = sqrt(100e-5);
bhat= 6.0776e6;

nr = 1e2;
ns = 1e4;
t0 = 0.1*ns;
y30 = zeros(nr,1);
y95 = zeros(nr,1);
tic;
for i=1:nr
	seed = i;
    [a,b]=challenger(ns+t0, 1, xdata, ydata, alphahat, alphasigma, betahat, betasigma, bhat, seed);
	eab30 = exp(a((t0+1):end)+b((t0+1):end)*30);
	eab95 = exp(a((t0+1):end)+b((t0+1):end)*95);
    y30(i) = mean(eab30./(1+eab30));
	y95(i) = mean(eab95./(1+eab95));
end
save('challenger.mat');