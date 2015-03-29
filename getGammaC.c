#include <mex.h>
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
    double *fx;
	double *cs;
	double *gamma, *gidx;
	int eta;
	int i, k, kk;
	int n, ng, nn;
    
    // Input
    n = mxGetM(prhs[0]);			// Number of steps
	fx = mxGetPr(prhs[0]);			// Function output
	gidx = mxGetPr(prhs[1]);		// Index of gamma terms to obtain
	ng = mxGetN(prhs[1]);			// Number of gamma terms
	eta = (int)mxGetPr(prhs[2])[0];	// Eta parameter
	nn = n - eta;
	
	// Output
	plhs[0] = mxCreateDoubleMatrix(1,ng,mxREAL);
	gamma = mxGetPr(plhs[0]);		// Output vector

	
	if(gidx[ng-1] > eta){
		mexErrMsgTxt("The maximal gamma index cannot exceed eta.");
		}
	
	//printf("%d,%d,%d,%d,%lf\n",n,eta,ng,nn,fx[0]);
	
	// Buuld cumulative sum vector
	cs = mxCalloc(nn+gidx[ng-1],sizeof(double));	
	cs[0] = fx[0];									
	for(i=1;i<nn+gidx[ng-1];i++){
		cs[i] = cs[i-1]+fx[i];
		}
	
	
	for(k=0;k<ng;k++){
		kk = gidx[k];
		gamma[k] = 0;
		for(i=0;i<nn;i++){
			gamma[k] += (fx[i]*fx[kk+i])/nn;
			}
		gamma[k] -= 0.5*(cs[nn-1]/nn)*(cs[nn-1]/nn);
		gamma[k] -= 0.5*((cs[nn+kk-1]-cs[kk-1])/nn)*((cs[nn+kk-1]-cs[kk-1])/nn);
		}
	
	mxFree(cs);
}
