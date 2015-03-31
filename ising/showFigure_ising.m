%% Compile all code
compile

%% Set basic parameters
addpath ..
N = 1e6;
t0 = ceil(0.1*N);
k = 10*ceil(N^(1/3));

%% CW, Glauber
latticeSize = 10; T = 2; h = 0;
seed = 1234;
[~,~,fx_cwglauber] = isingCW(1,latticeSize,N,t0,1,T,h,seed);
% Estimated parameters
gamma_cwglauber_est = getGamma(fx_cwglauber);
Vf_cwglauber_est = getVf(fx_cwglauber);
sigma_cwglauber_est = getSigmaNonasym(fx_cwglauber,k);
beta = 1/T;
C = latticeSize;
% Theoretical
tmix_cwglauber = 0.5*latticeSize*log((1-beta)^2*latticeSize)/(1-beta);
gamma_cwglauber = (1-beta)/latticeSize;
t0_cwglauber= floor(30*tmix_cwglauber);
load('./output/CW.mat','m','nSteps','nRelaxSteps');
NN = nSteps - nRelaxSteps;
% Simulation
[t,logp] = logTail(m,100);
% Chebyshev
logp_cheb = chebyshevtail(t,NN,sigma_cwglauber_est,Vf_cwglauber_est,gamma_cwglauber);
% Bernstein
logp_bern = bernsteintail(t,NN,sigma_cwglauber_est,Vf_cwglauber_est,gamma_cwglauber,C);
% Normal distribution
logp_norm = normasym(t,NN,sigma_cwglauber_est);


fh = logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm,true);
xlim([-0.4,0.4]);
%saveas(fh,'figures/cwglauber.eps');

fprintf(['Lattice size $%d$, $10^6$ runs, $N=10^5$, $t_0=%d$, $T=1/\\beta=%.2f$, $h=0$, $C=%d$, ', ...
	'$\\hat{\\sigma}^2 = %.5g$, $\\hat{V}_f = %.2f$, $\\hat{\\gamma}^{CW}_{Gl}=%0.2e$\n'],...
	latticeSize,t0_cwglauber,T,C,sigma_cwglauber_est,Vf_cwglauber_est,gamma_cwglauber);

%% CW, Metropolis
latticeSize = 10; T = 2; h = 0;
[~,~,fx_cwmetro] = isingCWMetropolis(1,latticeSize,N,t0,1,T,h,seed);
sigma_cwmetro_est = getSigmaNonasym(fx_cwmetro,k);
gamma_cwmetro_est = getGamma(fx_cwmetro);
Vf_cwmetro_est = getVf(fx_cwmetro);
tmix_cwmetro_est = 1/gamma_cwmetro_est;
t0_cwmetro = floor(30*tmix_cwmetro_est);
beta = 1/T;
C = latticeSize;
load('./output/CWMetropolis.mat','m','nSteps','nRelaxSteps');
NN = nSteps - nRelaxSteps;
% Simulation
[t,logp] = logTail(m,100);
% Chebyshev
logp_cheb = chebyshevtail(t,NN,sigma_cwmetro_est,Vf_cwmetro_est,gamma_cwmetro_est);
% Bernstein
logp_bern = bernsteintail(t,NN,sigma_cwmetro_est,Vf_cwmetro_est,gamma_cwmetro_est,C);
% Normal distribution
logp_norm = normasym(t,NN,sigma_cwmetro_est);

fh = logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm);
xlim([-0.3 0.3]);
%saveas(fh,'figures/cwmetro.eps');

fprintf(['Lattice size $%d$, $10^6$ runs, $N=10^5$, $t_0=%d$, $T=1/\\beta=%.2f$, $h=0$, $C=%d$, ', ...
	'$\\hat{\\sigma}^2 = %.5g$, $\\hat{V}_f = %.2f$, $\\hat{\\gamma}=%0.2e$\n'],...
	latticeSize,t0_cwmetro,T,C,sigma_cwmetro_est,Vf_cwmetro_est,gamma_cwmetro_est);

%% CW, Glauber, low temperature
latticeSize = 10; T = 0.5; h = 0;
p = 100; % Make a longer run here
[~,~,fx_cwlowtemp] = isingCW(1,latticeSize,N*p,t0*p,1,T,h,seed);
gamma_cwlowtemp_est = getGamma(fx_cwlowtemp);
sigma_cwlowtemp_est = getSigmaNonasym(fx_cwlowtemp,10*floor((N*p)^(1/3)));
Vf_cwlowtemp_est = getVf(fx_cwlowtemp);
tmix_cwlowtemp_est = 1/gamma_cwlowtemp_est;
t0_cwlowtemp = floor(30*tmix_cwlowtemp_est);
C = latticeSize;
load('./output/CWLowtemp.mat','m','nSteps','nRelaxSteps');
NN = nSteps - nRelaxSteps;
% Simulation
[t,logp] = logTail(m,100);

% Chebyshev
logp_cheb = chebyshevtail(t,NN,sigma_cwlowtemp_est,Vf_cwlowtemp_est,gamma_cwlowtemp_est);
% Bernstein
logp_bern = bernsteintail(t,NN,sigma_cwlowtemp_est,Vf_cwlowtemp_est,gamma_cwlowtemp_est,C);
% Normal distribution
logp_norm = normasym(t,NN,sigma_cwlowtemp_est);

fh = logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm);
%saveas(fh,'figures/cwlowtemp.eps');

fprintf(['Lattice size $%d$, $10^6$ runs, $N=10^5$, $t_0=%d$, $T=1/\\beta=%.2f$, $h=0$, $C=%d$, ', ...
	'$\\hat{\\sigma}^2 = %.5g$, $\\hat{V}_f = %.2f$, $\\hat{\\gamma}=%0.2e$\n'],...
	latticeSize,t0_cwlowtemp,T,C,sigma_cwlowtemp_est,Vf_cwlowtemp_est,gamma_cwlowtemp_est);

%% 1D, Glauber, random scan
latticeSize = 10; T = 2; h = 0;
[~,~,fx_1drand] = ising1D(1,latticeSize,N,t0,1,T,h,seed);
beta = 1/T;
% Theoretical
gamma_1drand = (2/latticeSize)*exp(-4*beta)/(1+exp(-4*beta));
tmix_1drand = (latticeSize/2)*(1+exp(4*beta))*log(4*latticeSize);
t0_1drand = floor(30*tmix_1drand);
% Estimated
gamma_1drand_est = getGamma(fx_1drand(:));
sigma_1drand_est = getSigmaNonasym(fx_1drand(:),k);
Vf_1drand_est = getVf(fx_1drand);

C = latticeSize;
load('./output/I1.mat','m','nSteps','nRelaxSteps');
NN = nSteps - nRelaxSteps;

% Simulation
[t,logp] = logTail(m,100);
% Chebyshev
logp_cheb = chebyshevtail(t,NN,sigma_1drand_est,Vf_1drand_est,gamma_1drand);
% Bernstein
logp_bern = bernsteintail(t,NN,sigma_1drand_est,Vf_1drand_est,gamma_1drand,C);
% Normal distribution
logp_norm = normasym(t,NN,sigma_1drand_est);

fh = logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm);
%saveas(fh,'figures/i1drand.eps');


fprintf(['Lattice size $%d$, $10^6$ runs, $N=10^5$, $t_0=%d$, $T=1/\\beta=%.2f$, $h=0$, $C=%d$, ', ...
	'$\\hat{\\sigma}^2 = %.5g$, $\\hat{V}_f = %.2f$, $\\hat{\\gamma}^{1D}_{Gl}=%0.2e$\n'],...
	latticeSize,t0_1drand,T,C,sigma_1drand_est,Vf_1drand_est,gamma_1drand);


%% 1D, Glauber, systematic scan
latticeSize = 10; T = 2; beta = 1/T; h = 0;
[~,~,fx_1dsyst] = ising1Dsyst(1,latticeSize,ceil(N/latticeSize),ceil(t0/latticeSize),1,T,h,seed);
% Theoretical
gamma_1dsyst = 8*exp(-4*beta)*(1+exp(-4*beta)) / ((1+3*exp(-4*beta))^2);
tmix_1dsyst = 0.25*(3+exp(4*beta))*log(4*latticeSize);
t0_1dsyst = floor(30*tmix_1dsyst);
% Estimated
gamma_1dsyst_est = getGamma(fx_1dsyst(:));
sigma_1dsyst_est = getSigmaNonasym(fx_1dsyst(:),ceil(k/latticeSize));
Vf_1dsyst_est = getVf(fx_1dsyst);

C = latticeSize;
load('./output/I1Syst.mat','m','nSteps','nRelaxSteps');
NN = nSteps - nRelaxSteps;
% Simulation
[t,logp] = logTail(m,100);

% Chebyshev
logp_cheb = chebyshev_nonrev(t,NN,sigma_1dsyst_est,Vf_1dsyst_est,gamma_1dsyst);
% Bernstein
logp_bern = bernstein_nonrev(t,NN,sigma_1dsyst_est,Vf_1dsyst_est,gamma_1dsyst,C);
% Normal distribution
logp_norm = normasym(t,NN,sigma_1dsyst_est);

fh = logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm);
xlim([-0.6 0.6]);
%saveas(fh,'figures/i1dsyst.eps');


fprintf(['Lattice size $%d$, $10^6$ runs, $N=10^5$, $t_0=%d$, $T=1/\\beta=%.2f$, $h=0$, $C=%d$, ', ...
	'$\\hat{\\sigma}^2 = %.5g$, $\\hat{V}_f = %.2f$, $(\\hat{\\gamma}_{\\mathrm{ps})^{1D}_{Gl.sys}=%0.2e$\n'],...
	latticeSize,t0_1dsyst,T,C,sigma_1dsyst_est,Vf_1dsyst_est,gamma_1dsyst);

%% CW, Glauber, sign of magnetization
latticeSize = 10; T = 2; h = 2;
[~,~,fx_cwmagsign] = isingCWMagsign(1,latticeSize,N,t0,1,T,h,seed);
% Estimated
gamma_cwmagsign_est = getGamma(fx_cwmagsign(:));
sigma_cwmagsign_est = getSigmaNonasym(fx_cwmagsign(:),k);
Vf_cwmagsign_est = getVf(fx_cwmagsign);
tmix_cwmagsign_est = 1/gamma_cwmagsign_est;
t0_cwmagsign = floor(30*tmix_cwmagsign_est);

beta = 1/T;
C = 1;
load('./output/CWMagsign.mat','m','nSteps','nRelaxSteps');
NN = nSteps - nRelaxSteps;
% Simulation
[t,logp] = logTail(m,100);

% Chebyshev
logp_cheb = chebyshevtail(t,NN,sigma_cwmagsign_est,Vf_cwmagsign_est,gamma_cwmagsign_est);
% Bernstein
logp_bern = bernsteintail(t,NN,sigma_cwmagsign_est,Vf_cwmagsign_est,gamma_cwmagsign_est,C);
% Normal distribution
logp_norm = normasym(t,NN,sigma_cwmagsign_est);

fh = logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm);
%saveas(fh,'figures/cwmagsign.eps');

fprintf(['Lattice size $%d$, $10^6$ runs, $N=10^5$, $t_0=%d$, $T=1/\\beta=%.2f$, $h=2$, $C=%d$, ', ...
	'$\\hat{\\sigma}^2 = %.5g$, $\\hat{V}_f = %.2f$, $\\hat{\\gamma}=%0.2e$\n'],...
	latticeSize,t0_cwmagsign,T,C,sigma_cwmagsign_est,Vf_cwmagsign_est,gamma_cwmagsign_est);
