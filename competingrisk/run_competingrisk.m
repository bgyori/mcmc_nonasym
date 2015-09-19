niters = 5000;
ns = 1e4;
t0 = 2e4;
probtla = zeros(1, niters);

mxdiv = 200;
mxt = 20;
dt = mxt/mxdiv;
t = linspace(mxt/mxdiv,mxt,mxdiv);
a = 1;

parfor j=1:niters
    [out,outend,data] = competingrisk_mex(ns+t0);
    [r,nstages] = size(out);
    for i=1:ns
        b1 = out(t0+i,nstages-3);
        sig1 = out(t0+i,nstages-2);
        b2 = out(t0+i,nstages-1);
        sig2 = out(t0+i,nstages);      
        
        h1 = b1*t.^(b1-1)./(sig1^b1*(1+(t/sig1).^b1));
        h2 = b2*t.^(b2-1)./(sig2^b2*(1+(t/sig2).^b2));
        s1 = 1./(1+(t/sig1).^b1);
        s2 = 1./(1+(t/sig2).^b2);
        
        pd1 = (h1+h2).*s1.*s2;
        pd0 = s1.*s2;
        spd1 = sum(pd1);
        spd0 = sum(pd0);
        psum = spd1+spd0;
        
        p = (pd1+pd0)/psum;
        
        Et=sum(p.*t);
        
        probtla(j) = probtla(j) + (Et>a);
    end
end
probtla = probtla/ns;
save('competingrisk_out.mat')