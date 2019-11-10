#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "NURBS.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Return the univariate NURBS basis functions and first derivatives to matlab
     * // All non-zero basis functions and derivatives at point [xi] are
     * // computed.
     * //
     * // We expect the function to be called as
     * // [R dRdxi]      = NURBS1DBasisDers(xi,p,knotU,weights)
     * //
     * //	xi           = where we want to interpolate
     * //	knotU        = knot vector
     * //	weights      = vector of weights
     * //
     * // Vinh Phu Nguyen, nvinhphu@gmail.com
     */
    
    if(nrhs != 4) mexErrMsgTxt("You fool! You haven't passed in 6 arguments to the function."
            "We expect it to be in the form [dRdxi dRdeta] = NURBSinterpolation(xi, p, q, knotU, knotV, weights)\n");
    
    /* First get the inputs */
    
    double *xi      = mxGetPr(prhs[0]);
    double *p_in    = mxGetPr(prhs[1]);
    int    p        = (int) *p_in;
    double *knotU   = mxGetPr(prhs[2]);
    int    numKnotU = mxGetN(prhs[2]);
    int    nU       = numKnotU - 1 - p - 1;
    int    noFuncs  = (p+1);
    double *weight  = mxGetPr(prhs[3]);
    int    numWeights = mxGetN(prhs[3]);
    
    double tol    = 100*DBL_EPSILON;
    
    if(fabs(xi[0]-knotU[numKnotU-1]) < tol)
        xi[0] = knotU[numKnotU-1] - tol;
    
    /* and evaluate the non-zero univariate B-spline basis functions
     * and first derivatives
     */
    
    double *N      = (double *)malloc(sizeof(double)*(p+1));
    double **dersN = init2DArray(nU+1, p+1);
    
    int    spanU   = FindSpan(nU, p, *xi, knotU);
    
    BasisFuns     (spanU, *xi, p, knotU, N);
    dersBasisFuns (spanU, *xi, p, nU, knotU, dersN);
    
    
    /* and create NURBS approximation */
    
    int i;
    
    /*
     * for(i=0;i<=p;i++){
     * printf("dNdxi= %f\n", dersN[1][i]);
     * printf("dNdet= %f\n", dersM[1][i]);
     * }*/
    
    int    uind  = spanU - p;
    double w     = 0.0;
    double dwdxi = 0.0;
    double wgt;
    
    for(i = 0; i <= p; i++)
    {       
        wgt    = weight[uind+i];        
        w     += N[i]        * wgt;
        dwdxi += dersN[1][i] * wgt;        
    }
        
    /* create output */
    
    plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,noFuncs,mxREAL);
    
    double *R      = mxGetPr(plhs[0]);
    double *dRdxi  = mxGetPr(plhs[1]);
    
    double fac;
    
    /*printf("uind= %d\n", uind);
     * printf("vind= %d\n", vind);
     * printf("nU+1= %d\n", nU+1);  */
    
    for(i = 0; i <= p; i++)
    {
        fac      = weight[uind+i]/(w*w);
        
        R[i]     = N[i]*fac*w;
        dRdxi[i] = (dersN[1][i]*w - N[i]*dwdxi) * fac;                       
    }
       
    free(N);
    free2Darray(dersN, (nU+1));
}







