#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "NURBS.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Return the 3D NURBS basis functions and first derivatives to matlab
     * // All non-zero basis functions and derivatives at point [xi,eta,zeta] are
     * // computed.
     * //
     * // We expect the function to be called as
     * // [R dRdxi dRdeta dRdzeta] = NURBS3DBasisDers(xi,p,q,r, knotU, ...
     * //                                               knotV,knotZ, weights)
     * //
     * //	xi                  = point, [xi eta zeta], where we want to interpolate
     * //	knotU, knotV, knotZ = knot vectors
     * //	weights             = vector of weights
     * //
     * // Vinh Phu Nguyen, nvinhphu@gmail.com
     */
    
    if(nrhs != 8) mexErrMsgTxt("You fool! You haven't passed in 6 arguments to the function."
            "We expect it to be in the form [dRdxi dRdeta] = NURBSinterpolation(xi,p,q,r,knotU,knotV,knotW, weights)\n");
    
    /* First get the inputs */
    
    double *xi      = mxGetPr(prhs[0]);
    double *p_in    = mxGetPr(prhs[1]);
    double *q_in    = mxGetPr(prhs[2]);
    double *r_in    = mxGetPr(prhs[3]);
    
    int    p        = (int) *p_in;
    int    q        = (int) *q_in;
    int    r        = (int) *r_in;
    
    double *knotU   = mxGetPr(prhs[4]);
    double *knotV   = mxGetPr(prhs[5]);
    double *knotW   = mxGetPr(prhs[6]);
    
    int    numKnotU = mxGetN(prhs[4]);
    int    numKnotV = mxGetN(prhs[5]);
    int    numKnotW = mxGetN(prhs[6]);
    
    int    nU       = numKnotU - 1 - p - 1;
    int    nV       = numKnotV - 1 - q - 1;
    int    nW       = numKnotW - 1 - r - 1;
    int    noFuncs  = (p+1)*(q+1)*(r+1);
    
    double *weight    = mxGetPr(prhs[7]);
    int    numWeights = mxGetN(prhs[7]);
    
    double tol    = 100*DBL_EPSILON;
    
    if(fabs(xi[0]-knotU[numKnotU-1]) < tol)
        xi[0] = knotU[numKnotU-1] - tol;
    
    if(fabs(xi[1]-knotV[numKnotV-1]) < tol)
        xi[1] = knotV[numKnotV-1] - tol;
    
    if(fabs(xi[2]-knotW[numKnotW-1]) < tol)
        xi[2] = knotW[numKnotW-1] - tol;
    
    /* and evaluate the non-zero univariate B-spline basis functions
     * and first derivatives
     */
    
    double *N      = (double *)malloc(sizeof(double)*(p+1));
    double *M      = (double *)malloc(sizeof(double)*(q+1));
    double *P      = (double *)malloc(sizeof(double)*(r+1));
    
    double **dersN = init2DArray(nU+1, p+1);
    double **dersM = init2DArray(nV+1, q+1);
    double **dersP = init2DArray(nW+1, r+1);
    
    int spanU = FindSpan(nU, p, xi[0], knotU);
    int spanV = FindSpan(nV, q, xi[1], knotV);
    int spanW = FindSpan(nW, r, xi[2], knotW);
    
    BasisFuns     (spanU, xi[0], p, knotU, N);
    BasisFuns     (spanV, xi[1], q, knotV, M);
    BasisFuns     (spanW, xi[2], r, knotW, P);
    
    dersBasisFuns (spanU, xi[0], p, nU, knotU, dersN);
    dersBasisFuns (spanV, xi[1], q, nV, knotV, dersM);
    dersBasisFuns (spanW, xi[2], r, nW, knotW, dersP);
    
    /* and create NURBS approximation */
    
    int i, j, k, c, kk;
    
    /*
     * for(i=0;i<=p;i++){
     * printf("dNdxi= %f\n", dersN[1][i]);
     * printf("dNdet= %f\n", dersM[1][i]);
     * }*/
    
    int uind = spanU - p;
    int vind, wind;
    
    double w     = 0.0;
    double dwdxi = 0.0;
    double dwdet = 0.0;
    double dwdze = 0.0;
    double wgt;
    
    for(k = 0; k <= r; k++)
    {
        wind = spanW - r + k;
        
        for(j = 0; j <= q; j++)
        {
            vind = spanV - q + j;
            
            for(i = 0; i <= p; i++)
            {
                c   = uind + i + (nU+1) * ( (nV+1)*wind + vind);
                wgt = weight[c];
                
                w     += N[i]        * M[j] * P[k] * wgt;
                dwdxi += dersN[1][i] * M[j] * P[k] * wgt;
                dwdet += dersM[1][j] * N[i] * P[k] * wgt;
                dwdze += dersP[1][k] * N[i] * M[j] * wgt;
            }
        }
    }
    
    /* create output */
    
    plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,noFuncs,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,noFuncs,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,noFuncs,mxREAL);
    
    double *R      = mxGetPr(plhs[0]);
    double *dRdxi  = mxGetPr(plhs[1]);
    double *dRdet  = mxGetPr(plhs[2]);
    double *dRdze  = mxGetPr(plhs[3]);
    
    uind = spanU - p;
    kk   = 0;
    
    double fac;
    double nmp;
    
    /*printf("uind= %d\n", uind);
     * printf("vind= %d\n", vind);
     * printf("nU+1= %d\n", nU+1);  */
    
    for(k = 0; k <= r; k++)
    {
        wind = spanW - r + k;
        
        for(j = 0; j <= q; j++)
        {
            vind = spanV - q + j;
            
            for(i = 0; i <= p; i++)
            {
                c        = uind + i + (nU+1) * ( (nV+1)*wind + vind);
                fac      = weight[c]/(w*w);
                nmp      = N[i]*M[j]*P[k];
                
                R[kk]     = nmp * fac * w;
                dRdxi[kk] = (dersN[1][i]*M[j]*P[k]*w - nmp*dwdxi) * fac;
                dRdet[kk] = (dersM[1][j]*N[i]*P[k]*w - nmp*dwdet) * fac;
                dRdze[kk] = (dersP[1][k]*N[i]*M[j]*w - nmp*dwdze) * fac;
                
                kk      += 1;
            }
        }
    }
    
    /*mexPrintf("\nWe have a knot vector with %d components\n", numWeights);
     * for(k = 0; k < numWeights; k++) mexPrintf("%2.2f\t", weight[k]);
     * mexPrintf("\n");*/
    
    free(N);
    free(M);
    free(P);
    free2Darray(dersN, (nU+1));
    free2Darray(dersM, (nV+1));
    free2Darray(dersP, (nW+1));
}







