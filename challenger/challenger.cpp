#include "mex.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

extern void _main();

double pi2 = atan(1.0)*8.0;

double randu(){
   return (double)(rand())/(double)(RAND_MAX);
}

double randn(){
  static int second;
  static double r1, r2;
  if(second == 0){
	r1 = randu(); r2 = randu();
	second = 1;
	return cos(pi2*r1)*sqrt(-2.0*log(r2));
	}
  else {
	second = 0;
	return sin(pi2*r1)*sqrt(-2.0*log(r2));
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
	if (nrhs < 9) mexErrMsgTxt("Requires 9 input arguments: ns, xdata, ydata, alphahat, alphasigma, betahat, betasigma, bhat, seed");
	
	int ns = (int)mxGetPr(prhs[0])[0];
	double *xdata = (double*)mxGetPr(prhs[1]);
	double *ydata = (double*)mxGetPr(prhs[2]);
	double alphahat = (double)mxGetPr(prhs[3])[0];
	double alphasigma = (double)mxGetPr(prhs[4])[0];
	double betahat = (double)mxGetPr(prhs[5])[0];
	double betasigma = (double)mxGetPr(prhs[6])[0];
	double bhat = (double)mxGetPr(prhs[7])[0];
	int seed = (int)mxGetPr(prhs[8])[0];
	
	int ndata = mxGetM(prhs[1])>mxGetN(prhs[1])?mxGetM(prhs[1]):mxGetN(prhs[1]);
	
	plhs[0] = mxCreateDoubleMatrix(ns,1,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(ns,1,mxREAL);
	double* xalpha = mxGetPr(plhs[0]);
	double* xbeta = mxGetPr(plhs[1]);
	
	srand(seed);
	
	mxArray *xa = mxCreateDoubleMatrix(ns,2,mxREAL);
	double* x = mxGetPr(xa);

	double xprop[2];
	double p,lp,lpold,r;
	

	// Initialize state at (ahat,bhat)
	x[0] = xalpha[0] = alphahat;
	x[ns] = xbeta[0] = betahat;
	lpold = lpri(alphahat,bhat) + llh(xdata,ydata,ndata,alphahat,betahat);
	// For each step
	for(int j=1;j<ns;j++){
		// Proposal
		xprop[0] = x[j-1]+alphasigma*randn();
		xprop[1] = x[ns+j-1]+betasigma*randn();
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
			}
		else {
			// Keep previous state
			x[j] = x[j-1];
			x[ns+j] = x[ns+j-1];
			}

		// Function to save
		xalpha[j] = x[j];
		xbeta[j] = x[ns+j];
		}
	mxDestroyArray(xa);
	return;
}
