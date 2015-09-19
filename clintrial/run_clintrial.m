compile;

t0 = 600;
nbiters = 1e6;
ns = 1000;
avarray = zeros(nbiters,1);
parfor it=1:nbiters
    f = clintrial_mex(ns+t0);
    avarray(it) = mean(f(t0+1:end));
end
save('clintrial_out.mat');