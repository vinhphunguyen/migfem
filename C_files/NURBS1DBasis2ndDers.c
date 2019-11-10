#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "NURBS.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Return the univariate NURBS basis functions and first/second derivatives to matlab
     *  All non-zero basis functions and derivatives at point [xi] are
     *  computed.
     * 
     *  We expect the function to be called as
     *  [R dRdxi dR2dxi]      = NURBS1DBasis2ndDers(xi,p,knotU,weights)
     * 
     * 	xi           = where we want to interpolate
     * 	knotU        = knot vector
     * 	weights      = vector of weights
     * 
     *  Vinh Phu Nguyen, nvinhphu@gmail.com
     */
    
    if(nrhs != 4) mexErrMsgTxt("You fool! You haven't passed in 4 arguments to the function."
            "We expect it to be in the form [R dR2dxi] = NURBSinterpolation(xi, p, q, knotU, knotV, weights)\n");
    
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
    
    /*printf("vind = %d\n", numWeights);printf("vind = %f\n", xi[0]);*/
    
    /* and create NURBS approximation */
    
    int i;
    
    /*
     * for(i=0;i<=p;i++){
     * printf("dNdxi= %f\n", dersN[1][i]);
     * printf("dNdet= %f\n", dersM[1][i]);
     * }*/
    
    int    uind   = spanU - p;
    double w      = 0.0; /* N_iw_i */
    double dwdxi  = 0.0; /* first derivative of w */
    double d2wdxi = 0.0; /* second derivative of w */
    double wi;           /* wi */
    
    for(i = 0; i <= p; i++)
    {       
        wi      = weight[uind+i];        
        w      += N[i]        * wi;
        dwdxi  += dersN[1][i] * wi;        
        d2wdxi += dersN[2][i] * wi;        
    }
        
    /* create output */
    
    plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,noFuncs,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,noFuncs,mxREAL);
    
    double *R      = mxGetPr(plhs[0]);
    double *dRdxi  = mxGetPr(plhs[1]);
    double *dR2dxi = mxGetPr(plhs[2]);
    
    double fac;
    
    /*printf("uind= %d\n", uind);
     * printf("vind= %d\n", vind);
     * printf("nU+1= %d\n", nU+1);  */
    
    for(i = 0; i <= p; i++)
    {
        wi        = weight[uind+i]; 
        fac       = wi/(w*w);
        R[i]      = N[i]*wi/w;
        dRdxi[i]  = (dersN[1][i]*w - N[i]*dwdxi) * fac;  
        dR2dxi[i] = wi*(dersN[2][i]/w - 2*dersN[1][i]*dwdxi/w/w - N[i]*d2wdxi/w/w + 2*N[i]*dwdxi*dwdxi/w/w/w) ;                       
    }
       
    free(N);
    free2Darray(dersN, (nU+1));
}







