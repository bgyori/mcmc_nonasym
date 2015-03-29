
#include "mex.h"
#include "mexFunc.h"
#include "ising1DMain.h"


using namespace std;

extern void _main();

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
	double *magMean, *magVar;
	std::vector<double> mags;

	int nRuns, latticeSize, nSteps, nRelaxSteps, nStepInterval,getAllMag;
	double T, H;
	int seed;

	if (nrhs < 8){
		mexErrMsgTxt("Requires 8 input arguments.");
	}

	nRuns = (int)mxGetPr(prhs[0])[0];
	latticeSize = (int)mxGetPr(prhs[1])[0];
	nSteps = (int)mxGetPr(prhs[2])[0];
	nRelaxSteps = (int)mxGetPr(prhs[3])[0];
	nStepInterval = (int)mxGetPr(prhs[4])[0];
	T = (double)mxGetPr(prhs[5])[0];
	H = (double)mxGetPr(prhs[6])[0];
	seed = (int)mxGetPr(prhs[7])[0];
	if (nlhs == 3){
		getAllMag = 1;
		printf("Get all mag: %d\n",getAllMag);
	} else {
		getAllMag = 0;
	}
	
	if(nSteps<0) {
		mexErrMsgTxt("Negative step size.");
		}

	plhs[0] = mxCreateDoubleMatrix(nRuns,1,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(nRuns,1,mxREAL);
	magMean = mxGetPr(plhs[0]);
	magVar = mxGetPr(plhs[1]);


	if (getAllMag){
		magIsing1D(1,latticeSize,nSteps,nRelaxSteps,nStepInterval,T,H,seed,0,0,mags);
		plhs[2] = getMexArray(mags);
	} else {
		magIsing1D(nRuns,latticeSize,nSteps,nRelaxSteps,nStepInterval,T,H,seed,magMean,magVar,mags);
	}

	return;
}
