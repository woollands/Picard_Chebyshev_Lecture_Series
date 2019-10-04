/*=================================================================
 * EGMGravMex.c - example used to illustrate how to fill an mxArray
 *
 * Calculates gravity according to the EGM 2008 Spherical 
 * Harmonic Model
 *
 * Input:   Position Vector, Degree
 * Output:  Gravity Vector
 *
 * Copyright (c) TAMU LASR Lab, Austin Probe 2014
 *	
 *=================================================================*/
#include "mex.h"
#include "./EGM2008.c"


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
        //mexPrintf("\nGot Here.\n");
    
    double *inMatrix;               /* 1xN input matrix */
    int deg;                        /* Input degree */
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }

    
    /* make sure the first input argument is type double */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    
    /* check that number of rows in first input argument is 1 */
    if(mxGetM(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }

    if (!mxIsNumeric(prhs[1])) {
        mexErrMsgIdAndTxt( "MATLAB:timestwoalt:invalidInputType",
                "Argument 2 must be numeric.");
      } else if (mxGetNumberOfElements(prhs[1]) != 1 || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt( "MATLAB:timestwoalt:inputNotRealScalar",
                "Argument 2 must be non-complex scalar.");
      }
    
    /* create a pointer to the real data in the input matrix  */
    inMatrix = mxGetPr(prhs[0]);

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[0]);

    /* get the input degree */
    deg = mxGetScalar(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* call EGM Function */
    EGM2008(inMatrix, outMatrix, deg);
    
}
