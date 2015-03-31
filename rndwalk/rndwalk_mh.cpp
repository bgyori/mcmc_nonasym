#include "mex.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "rng_mersenne.h"

extern void _main();


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
		
	rng_mersenne rng(seed);
	std::cout << "Rand test: " << rng.randu() << std::endl;
	
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
			proposal[k] = rng.randn();
			xprop[k] = x[k*ns+j-1]+sigmap[k]*proposal[k];
			llh -= xprop[k]*xprop[k]/(2.0*sigmaf*sigmaf);
			}
		// Diff likelihood
		dllh = llh-llhold;
		// If new position is better or worse but randomly accept
		if((dllh>=0) || (exp(dllh)>rng.randu())){
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
