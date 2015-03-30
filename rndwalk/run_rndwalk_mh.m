compile;

ns = 1e4; % Number of samples
nr = 1e6; % Number of runs
sigmaf = 1; % Sigma of the function

%matlabpool open 4;

% Proposal sigmas
sigmaps = {1,0.5,[1,1,1,1,1],[0.7,0.5,0.4,0.6,1]};

% For each version of the proposal sigma
for i=1:length(sigmaps)
	sigmap = sigmaps{i};
	d = length(sigmap); % Number of dimensions
	m = zeros(nr,1); % Mean for each run
	parfor j=1:nr
		seed = j;
		% Run the random walk
        f = rndwalk_mh(2*ns,1,d,sigmap,sigmaf,3.0,seed);
        % Store the mean for each run, discarding initial ns samples
		m(j) = mean(f(ns+1:end));
	end
	fname = sprintf('output/rndwalk_mh%d.mat',i);
	save(fname,'m','ns','sigmaf','sigmap');
end

%matlabpool close;