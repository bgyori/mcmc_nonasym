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

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
	if (nrhs < 4) mexErrMsgTxt("Requires 6 input arguments: ns, sigmap, sigmaf, cutoff, seed");
	
	int ns = (int)mxGetPr(prhs[0])[0];
	double *sigmap = (double*)mxGetPr(prhs[1]);
	double sigmaf = (double)mxGetPr(prhs[2])[0];
	double cutoff = (double)mxGetPr(prhs[3])[0];
	int seed = (int)mxGetPr(prhs[4])[0];
	
	int d = mxGetN(prhs[1]);
	
	plhs[0] = mxCreateDoubleMatrix(ns,1,mxREAL);
	double* f = mxGetPr(plhs[0]);
		
	srand(seed);
	std::cout << "Rand test: " << rand() << std::endl;
	
	mxArray *xa = mxCreateDoubleMatrix(ns,d,mxREAL);
	double* x = mxGetPr(xa);
	if(nlhs > 1){
		plhs[1] = xa;
		}

	double proposal[d];
	double xprop[d];
	double p,dllh,llh,llhold=0.0;
	// Initialize state at (0,0,...,0)
	for(int k=0;k<d;k++){
		x[k*ns] = 0.0;
		}
	// For each step
	for(int j=1;j<ns;j++){
		llh = 0.0;
		for(int k=0;k<d;k++){
			proposal[k] = randn();
			xprop[k] = x[k*ns+j-1]+sigmap[k]*proposal[k];
			llh -= xprop[k]*xprop[k]/(2.0*sigmaf*sigmaf);
			}
		// Diff likelihood
		dllh = llh-llhold;
		// If new position is better or worse but randomly accept
		if((dllh>=0) || (exp(dllh)>randu())){
			// Accept proposal
			for(int k=0;k<d;k++){
				x[k*ns+j] = xprop[k];
				}
			llhold = llh;
			//std::cout << "accept" << std::endl;
			}

		else {
			// Keep previous state
			for(int k=0;k<d;k++){
				x[k*ns+j] = x[k*ns+j-1];
				}
			//std::cout << "reject" << std::endl;
			}

		// Function to save
		f[j] = (x[j]>cutoff);
		}
	if(nlhs<2){
		mxDestroyArray(xa);
		}
	return;
}
