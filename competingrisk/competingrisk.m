function [out, outend, data] = competingrisk(ns)
    assert(isa(ns, 'double'));
    
	% Load existing data for reproducibility
    D = coder.load('data.mat');
    data = D.data;
	
    nstages = 5;
    pararray = rand(nstages,4)*3.6 + 1;
    out = zeros(ns,4*nstages);
    outend = zeros(ns,4);
    acceptratear = zeros(nstages,2);
    for i=1:ns
        [pararray, acc] = mcmcstep(pararray, nstages, data);
        acceptratear = acceptratear + acc;
        for j = 1:nstages
            out(i,(j-1)*4+1:j*4) = pararray(j,:);
        end
        outend(i,1:4) = pararray(nstages,:);
	end
end

% Log-posterior evaluation
function logp = mu(par, data, betak)
	b1 = par(1);
	sig1 = par(2);
	b2 = par(3);
	sig2 = par(4);

	if (b1<=0 || b2<=0 || sig1<=0 || sig2<=0 || ...
			b1>5 || b2>5 || sig1>5 || sig2>5)
		logp=-inf;
		return;
	end

	t = data(:,1);
	d = data(:,1);
	h1 = b1*t.^(b1-1)./(sig1^b1*(1+(t/sig1).^b1));
	h2 = b2*t.^(b2-1)./(sig2^b2*(1+(t/sig2).^b2));
	s1 = 1./(1+(t/sig1).^b1);
	s2 = 1./(1+(t/sig2).^b2);
	logp = sum(d.*log(h1 + h2) + log(s1) + log(s2));
	logp = betak*logp;

end

% One step of the Markov chain
function [newpar, acceptrate] = mcmcstep(pararray,nstages,data)
	ms = 1 * fliplr(linspace(0.2,1,nstages))' * ones(1,4);
	newpar = pararray;
	acceptrate = zeros(nstages,2);

	betaar = linspace(0,1,nstages);

	k = floor(rand*(nstages-1)) + 1;
	betak = betaar(k);
	betakp = betaar(k+1);
	q = mu(newpar(k+1,:), data, betak) + ...
		mu(newpar(k,:), data, betakp) - ...
		(mu(newpar(k,:), data, betak) + ...
		mu(newpar(k+1,:), data, betakp));

	acceptrate(k,1) = acceptrate(k,1) + 1;
	r = rand();
	if log(r) < q
		acceptrate(k,2) = acceptrate(k,2) + 1;
		z = newpar(k,:);
		newpar(k,:) = newpar(k+1,:);
		newpar(k+1,:) = z;
	end

	for k = 1:nstages
		betak = betaar(k);
		newp = pararray(k,:);
		newp(1) = newp(1) + randn*ms(k,1);
		newp(2) = newp(2) + randn*ms(k,2);
		newp(3) = newp(3) + randn*ms(k,3);
		newp(4) = newp(4) + randn*ms(k,4);
		lograt = mu(newp,data,betak) - ...
				mu(pararray(k,:),data,betak);
		r = rand();
		if log(r) < lograt
			newpar(k,:) = newp;    
		end
	end

	k = floor(rand*(nstages-1)) + 1;
	betak = betaar(k);
	betakp = betaar(k+1);
	q = mu(newpar(k+1,:), data, betak) + ...
		mu(newpar(k,:), data, betakp) - ...
		(mu(newpar(k,:), data, betak) + ...
		mu(newpar(k+1,:), data, betakp));

	acceptrate(k,1) = acceptrate(k,1) + 1;
	r = rand;
	if log(r) < q
		acceptrate(k,2) = acceptrate(k,2) + 1;
		z = newpar(k,:);
		newpar(k,:) = newpar(k+1,:);
		newpar(k+1,:) = z;
	end
end

% Generate new data
function data = gendata(ns)
	b1 = 1;
	sig1 = 1;
	b2 = 4;
	sig2 = 4;

	mxdiv = 10000;
	mxt = 4.2;
	dt = mxt/mxdiv;
	t = linspace(0,mxt,mxdiv);

	h1 = b1*t.^(b1-1)./(sig1^b1*(1+(t/sig1).^b1));
	h2 = b2*t.^(b2-1)./(sig2^b2*(1+(t/sig2).^b2));
	s1 = 1./(1+(t/sig1).^b1);
	s2 = 1./(1+(t/sig2).^b2);

	pd1 = (h1+h2).*s1.*s2;
	pd0 = s1.*s2;
	psum = sum(pd1+pd0);
	probd1 = sum(pd1)/psum;
	probd0 = sum(pd0)/psum;

	pd1 = abs(pd1/sum(pd1));
	pd0 = abs(pd0/sum(pd0));

	cdf1 = cumsum(pd1);
	M = length(cdf1);
	xx = linspace(0,1,M);
	invcdf1 = interp1(cdf1,linspace(0,mxt,M),xx);


	cdf0 = cumsum(pd0);
	M = length(cdf0);
	xx = linspace(0,1,M);
	invcdf0 = interp1(cdf0,linspace(0,mxt,M),xx);

	data = zeros(ns,2);
	for i = 1:ns
		r = rand;
		if r < probd1
			data(i,1) = invcdf1(floor(rand*M)+1);
			data(i,2) = 1;
		else
			data(i,1) = invcdf0(floor(rand*M)+1);
			data(i,2) = 0;
		end
	end
end