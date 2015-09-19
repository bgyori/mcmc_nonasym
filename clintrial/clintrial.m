function [f, allvariables] = clintrial(ns)
    assert(isa(ns,'double'));
	coder.extrinsic('betarnd')
	coder.extrinsic('betacdf')
	
    % Data
    % MYX = [X=1/Y=1, X=0/Y=1; X=1/Y=0, X=0/Y=0]
    MYX = [10, 1; 5, 4];
    % MY = [Y=0, Y=1]
    MY = [10, 10];
    nx = sum(sum(MYX));
    ny = sum(MY);
    nz = 10;
    nmax = 100;

    % Hyperparameters
    % gamma0 ~ Beta(a0,b0)
    a0 = 1;
    b0 = 1;
    % gamma_1 ~ Beta(a1,b1)
    a1 = 1;
    b1 = 1;
    % P(X=1) = p
    % p ~ Beta(ax, bx)
    ax = 1;
    bx = 1;

    %nburn = 1e3;
    %nsamp = 1e4;
    gamma0 = zeros(ns+1,1);
    gamma1 = zeros(ns+1,1);
    xvec = zeros(ns+1,nx+ny+nz);
    p = zeros(ns+1,1);
    succ = zeros(ns,1);
    succ_pred = zeros(ns,1);
    gamma0(1) = rand*(1-2/168)+1/168;
    gamma1(1) = rand*(1-2/168)+1/168;
    p(1) = rand*(1-2/408)+1/408;
    xvec(1,:)=rand(1,nx+ny+nz)<1/2;
    for i=1:ns
        % Generate gamma0
        gamma0(i+1) = betarnd(a0 + MYX(2,1), b0 + MYX(2,2));
        % Generate gamma1
        gamma1(i+1) = betarnd(a1 + MYX(1,1), b1 + MYX(1,2));
        % Fixed values of X
        xvec(i+1, 1:nx) = [ones(1, sum(MYX(:,1))), ...
                            zeros(1, sum(MYX(:,2)))];
        % Predict X based on observed Y
        xvec(i+1, (nx+1):(nx+MY(1))) = ...
                    (rand(1, MY(1)) < gamma0(i+1));
        xvec(i+1, (nx+MY(1)+1):(nx+ny)) = ...
                    (rand(1, MY(2)) < gamma1(i+1));
        % Predict X based on no data
        xvec(i+1, (nx+ny+1):(nx+ny+nz)) = ...
                    (rand(1, nz) < p(i));
        nxpos = sum(xvec(i+1,:)==1);
        nxneg = sum(xvec(i+1,:)==0);
        % Sample p from full conditional
        p(i+1) = betarnd(ax + nxpos, bx + nxneg);
        
        if nargout > 1
            % Predict X based on no data
            xvec(i+1, (nx+ny+1):(nx+ny+nz)) = ...
                        (rand(1, nz) < p(i));
            nxpos = sum(xvec(i+1,:)==1);
            nxneg = sum(xvec(i+1,:)==0);    
            % Predict X based on observed Y
            xvec(i+1, (nx+1):(nx+MY(1))) = ...
                        (rand(1, MY(1)) < gamma0(i+1));
            xvec(i+1, (nx+MY(1)+1):(nx+ny)) = ...
                        (rand(1, MY(2)) < gamma1(i+1));
            % Generate gamma1
            gamma1(i+1) = betarnd(a1 + MYX(1,1), b1 + MYX(1,2));
            % Generate gamma0
            gamma0(i+1) = betarnd(a0 + MYX(2,1), b0 + MYX(2,2));
		end
        % Success probability at current sample size
		tmp = 0;
		tmp = betacdf(0.5, nxpos, nxneg);
        succ(i) = 1-tmp;
        % Success probability at maximal sample size
        pos_pred = sum(rand(1,nx+ny+nz) < p(i+1));
        nxpos_pred = nxpos + pos_pred;
        nxneg_pred = nxneg + (nx+ny+nz-pos_pred);
		tmp = betacdf(0.5, ax + nxpos_pred, ax + nxneg_pred);
        succ_pred(i) = 1-tmp;
    end
    
    f = (succ_pred > 0.5);
    
    if nargout > 1
        allvariables=[xvec gamma0 gamma1 p];
    end
end