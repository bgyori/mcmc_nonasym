function abtomb=challenger(mxiter)
assert(isa(mxiter,'double'));

x=[53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81];
y=[1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0];
gm=0.577216;
% banana=@(ab)-logL(ab(1),ab(2),x,y);
% [ab,minab]=fminsearch(banana,[10,0.2])
% alphahat=ab(1)
% betahat=ab(2)
% eps=1e-7;
% res=(banana([alphahat,betahat])-banana([alphahat,betahat+eps]))-(banana([alphahat,betahat-eps])-banana([alphahat,betahat]));
% res=res/eps^2;
% sigmabsq=-1/res;
% res=(banana([alphahat+eps,betahat])-banana([alphahat,betahat]))-(banana([alphahat,betahat])-banana([alphahat-eps,betahat]));
% res=res/eps^2;
% sigmaasq=1/res;
% 
% bhat=exp(alphahat+gm);

alphahat=15.0429;
betahat= -0.2322;
%sigmabsq=6.5952e-05;
sigmabsq=100e-5;
%sigmaasq=0.3311;
sigmaasq=4;
bhat= 6.0776e+06;

atomb=zeros(1,mxiter);
btomb=zeros(1,mxiter);
atomb(1)=alphahat;
btomb(1)=betahat;


for(iter=1:mxiter-1)
    [alpha,beta]=ggen(atomb(iter),btomb(iter),sigmaasq,sigmabsq);
    r=rand;
    %rat=exp(logL(alpha,beta,x,y)-logL(atomb(iter),btomb(iter),x,y)+logg(atomb(iter),btomb(iter),betahat,bhat,sigmabsq)-logg(alpha,beta,betahat,bhat,sigmabsq));
    rat=exp(logprob(alpha,beta,bhat,x,y)-logprob(atomb(iter),btomb(iter),bhat,x,y));
    if(r<rat) atomb(iter+1)=alpha;btomb(iter+1)=beta;
    else
        atomb(iter+1)=atomb(iter);btomb(iter+1)=btomb(iter);
    end
end

abtomb=[atomb;btomb]';
end

function p=L(a,b,x,y)
n=23;
p=1;
for(i=1:n)
    if(y(i)==1)
    p=p*(exp(a+b*x(i))/(1+exp(a+b*x(i))));
    else
    p=p*(1/(1+exp(a+b*x(i))));
    end
end
end

function lp=logL(a,b,x,y)
n=23;
lp=0;
for(i=1:n)
    if(y(i)==1)
    lp=lp+log((exp(a+b*x(i))/(1+exp(a+b*x(i)))));
    else
    lp=lp+log(1/(1+exp(a+b*x(i))));
    end
end
end




function p=prior(a,bhat)
    p=(1/bhat)*exp(a)*exp(-exp(a)/bhat);
end

function lp=logprior(a,bhat)
    lp=log(1/bhat)+a-exp(a)/bhat;
end

function p=prob(a,b,bhat,x,y)
p=L(a,b,x,y)*prior(a,bhat);
end

function lp=logprob(a,b,bhat,x,y)
    lp=logprior(a,bhat)+logL(a,b,x,y);
end


function p=g(a,b,betahat,bhat,sigmabsq)
p=exp(-(b-betahat)^2/(2*sigmabsq));
end

function p=logg(a,b,betahat,bhat,sigmabsq)
p=-(b-betahat)^2/(2*sigmabsq);
end

function [alpha,beta]=ggen(a,b,sigmaasq,sigmabsq)
alpha=a+randn*sqrt(sigmaasq);
beta=b+randn*sqrt(sigmabsq);
end

% function abarray=challenger(mxiter)
% assert(isa(mxiter,'double'));
% 
% x=[53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81];
% y=[1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0];
% gm=0.577216;
% alphahat=15.0429;
% betahat= -0.2322;
% sigmabsq=100e-5;
% sigmaasq=4;
% bhat= 6.0776e+06;
% 
% aarray=zeros(1,mxiter);
% barray=zeros(1,mxiter);
% aarray(1)=alphahat;
% barray(1)=betahat;
% 
% 
% for(iter=1:mxiter-1)
%     [alpha,beta]=ggen(aarray(iter),barray(iter),sigmaasq,sigmabsq);
%     r=rand;
%     rat=exp(logprob(alpha,beta,bhat,x,y)-logprob(aarray(iter),barray(iter),bhat,x,y));
%     if(r<rat) aarray(iter+1)=alpha; barray(iter+1)=beta;
%     else
%         aarray(iter+1)=aarray(iter); barray(iter+1)=barray(iter);
%     end
% end
% 
% abarray=[aarray;barray]';
% end
% 
% function lp=logL(a,b,x,y)
% n=23;
% lp=0;
% for(i=1:n)
%     if(y(i)==1)
%     lp=lp+log((exp(a+b*x(i))/(1+exp(a+b*x(i)))));
%     else
%     lp=lp+log(1/(1+exp(a+b*x(i))));
%     end
% end
% end
% 
% 
% function lp=logprior(a,bhat)
%     lp=log(1/bhat)+a-exp(a)/bhat;
% end
% 
% function lp=logprob(a,b,bhat,x,y)
%     lp=logprior(a,bhat)+logL(a,b,x,y);
% end
% 
% function p=logg(a,b,betahat,bhat,sigmabsq)
% p=-(b-betahat)^2/(2*sigmabsq);
% end
% 
% function [alpha,beta]=ggen(a,b,sigmaasq,sigmabsq)
% alpha=a+randn*sqrt(sigmaasq);
% beta=b+randn*sqrt(sigmabsq);
% end

% function [alpha, beta] = challenger(ns)
% 	assert(isa(ns,'double'));
% 
% 	x = [53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81];
% 	y = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0];
% 	gm = 0.577216;
% 	alphahat = 15.0429;
% 	betahat = -0.2322;
% 	sigmabsq = 1e-3;
% 	sigmaasq = 4;
% 	bhat = 6.0776e+06;
% 
% 	alpha = zeros(1,ns);
% 	beta = zeros(1,ns);
% 	alpha(1) = alphahat;
% 	beta(1) = betahat;
% 
% 
% 	for i=1:ns-1
% 		[a,b] = ggen(alpha(i),beta(i),sigmaasq,sigmabsq);
% 		r = rand;
% 		rat = exp(logprob(a,b,bhat,x,y) - ...
% 				logprob(alpha(i),beta(i),bhat,x,y));
% 		if(r<rat) 
% 			alpha(i+1) = a; 
% 			beta(i+1) = b;
% 		else
% 			alpha(i+1) = alpha(i); 
% 			beta(i+1) = beta(i);
% 		end
% 	end
% end
% 
% function lp = logL(a,b,x,y)
% 	n = 23;
% 	lp = 0;
% 	for i=1:n
% 		if y(i)==1
% 			lp = lp + log((exp(a+b*x(i))/(1+exp(a+b*x(i)))));
% 		else
% 			lp = lp + log(1/(1+exp(a+b*x(i))));
% 		end
% 	end
% end
% 
% function lp = logprior(a,bhat)
% 	lp = log(1/bhat) + a - exp(a)/bhat;
% end
% 
% function lp = logprob(a,b,bhat,x,y)
% 	lp = logprior(a,bhat) + logL(a,b,x,y);
% end
% 
% function p = logg(a,b,betahat,bhat,sigmabsq)
% 	p = -(b-betahat)^2/(2*sigmabsq);
% end
% 
% function [alpha,beta] = ggen(a,b,sigmaasq,sigmabsq)
% 	alpha = a + randn*sqrt(sigmaasq);
% 	beta = b + randn*sqrt(sigmabsq);
% end