#include "mex.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

extern void _main();

double randu(){
   return (double)(rand())/(double)(RAND_MAX);
}

void randn(double* r, int n){
  double r1, r2;
  double pi2 = atan(1.0)*8.0;
  for(int i=0;i<n;i+=2){
	r1 = randu(); r2 = randu();
	r[i] = cos(pi2*r1)*sqrt(-2.0*log(r2)); 
	if(i+1>=n)
		break;
	r[i+1] = sin(pi2*r1)*sqrt(-2.0*log(r2));
	}
}

double llh(double *xdata,double *ydata,int ndata,double alpha, double beta){
	double llh = 0.0;
	for(int i=0;i<ndata;i++){
		double eab = exp(alpha+beta*xdata[i]);
		if (ydata[i]==0)
			llh += -log(1+eab);
		else 
			llh += log(eab)-log(1+eab);
		}
	return llh;
	}

double lpri(double alpha, double b){
	return log(1/b) + alpha - exp(alpha)/b;
	}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
	if (nrhs < 10) mexErrMsgTxt("Requires 10 input arguments: ns, nr, xdata, ydata, alphahat, alphasigma, betahat, betasigma, bhat, seed");
	
	int ns = (int)mxGetPr(prhs[0])[0];
	int nr = (int)mxGetPr(prhs[1])[0];
	double *xdata = (double*)mxGetPr(prhs[2]);
	double *ydata = (double*)mxGetPr(prhs[3]);
	double alphahat = (double)mxGetPr(prhs[4])[0];
	double alphasigma = (double)mxGetPr(prhs[5])[0];
	double betahat = (double)mxGetPr(prhs[6])[0];
	double betasigma = (double)mxGetPr(prhs[7])[0];
	double bhat = (double)mxGetPr(prhs[8])[0];
	int seed = (int)mxGetPr(prhs[9])[0];
	
	int ndata = mxGetM(prhs[2])>mxGetN(prhs[2])?mxGetM(prhs[2]):mxGetN(prhs[2]);
	int d = 2;
	
	plhs[0] = mxCreateDoubleMatrix(ns,nr,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(ns,nr,mxREAL);
	double* xalpha = mxGetPr(plhs[0]);
	double* xbeta = mxGetPr(plhs[1]);
	
	srand(seed);
	//std::cout << "Rand test: " << rand() << std::endl;
	
	mxArray *xa = mxCreateDoubleMatrix(ns,d,mxREAL);
	double* x = mxGetPr(xa);

	double proposal[d];
	double xprop[d];
	double p,lp,lpold,r;
	
	// For each run
	for(int i=0;i<nr;i++){
		// Initialize state at (ahat,bhat)
		x[0] = xalpha[0] = alphahat;
		x[ns] = xbeta[0] = betahat;
		lpold = lpri(alphahat,bhat) + llh(xdata,ydata,ndata,alphahat,betahat);
		// For each step
		for(int j=1;j<ns;j++){
			// Proposal
			randn(proposal,d);
			xprop[0] = x[j-1]+alphasigma*proposal[0];
			xprop[1] = x[ns+j-1]+betasigma*proposal[1];
			// Log posterior
			lp = lpri(xprop[0],bhat) + llh(xdata,ydata,ndata,xprop[0],xprop[1]);
			// Log posterior ratio
			r = lp -lpold;
			
			// If new position is better or worse but randomly accept
			if((r>=0) || (exp(r)>randu())){
				// Accept proposal
				x[j] = xprop[0];
				x[ns+j] = xprop[1];
				lpold = lp;
				//std::cout << "accept" << std::endl;
				}
			else {
				// Keep previous state
				x[j] = x[j-1];
				x[ns+j] = x[ns+j-1];
				//std::cout << "reject" << std::endl;
				}
			
			// Function to save
			xalpha[i*ns+j] = x[j];
			xbeta[i*ns+j] = x[ns+j];
			}
		}
	mxDestroyArray(xa);
	return;
}
