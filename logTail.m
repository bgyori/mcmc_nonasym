function [t,logp] = logTail(m,N)
	if nargin < 2
		N = 100;
	end
	m = m(:)';
	nSample = length(m);
	meanM = mean(m);
	mShift = m - meanM;
	
	rangePos = max(m)-meanM;
	rangeNeg = meanM-min(m);
	ntPos = floor(rangePos/(rangePos+rangeNeg)*N);
	ntNeg = N - ntPos;
	
	tPos = rangePos*(0:ntPos)/ntPos;
	tNeg = -rangeNeg*(ntNeg:-1:0)/ntNeg;
	
	if(nSample <= 1e6)
		tailPos = (log(sum(repmat(mShift,ntPos+1,1)>repmat(tPos',1,nSample),2)/nSample))';
		tailNeg = (log(sum(repmat(mShift,ntNeg+1,1)<=repmat(tNeg',1,nSample),2)/nSample))';
	else
		for i=1:ntPos+1
			tailPos(i) = log(sum(mShift>tPos(i))/nSample)';
		end
		for i=1:ntNeg+1
			tailNeg(i) = log(sum(mShift<=tNeg(i))/nSample)';
		end
	end
	
	logp = [tailNeg,tailPos];
	t = [tNeg,tPos];
end