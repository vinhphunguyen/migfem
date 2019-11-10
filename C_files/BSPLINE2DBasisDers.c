#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "NURBS.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Return the NURBS basis functions and first derivatives to matlab
	 All non-zero basis functions and derivatives at point [xi,eta] are
     computed.

	 We expect the function to be called as 
	 [R dRdxi dRdeta] = NURBSinterpolation(xi, p, q, knotU, ...
                            knotV, weights)
	
		xi           = point, [xi eta], where we want to interpolate
		knotU, knotV = knot vectors
		weights      = vector of weights 
    
     Vinh Phu Nguyen, nvinhphu@gmail.com
    */
	
	if(nrhs != 5) mexErrMsgTxt("You fool! You haven't passed in 6 arguments to the function."
		"We expect it to be in the form [dRdxi dRdeta] = NURBSinterpolation(xi, p, q, knotU, knotV)\n");
			
	/* First get the inputs */
	
	double *xi      = mxGetPr(prhs[0]);	
	double *p_in    = mxGetPr(prhs[1]);
    double *q_in    = mxGetPr(prhs[2]);
    
	int    p        = (int) *p_in;
	int    q        = (int) *q_in;
    
	double *knotU   = mxGetPr(prhs[3]);
    double *knotV   = mxGetPr(prhs[4]);
    
	int    numKnotU = mxGetN(prhs[3]);
    int    numKnotV = mxGetN(prhs[4]);
    
	int    nU       = numKnotU - 1 - p - 1;
    int    nV       = numKnotV - 1 - q - 1;
    int    noFuncs  = (p+1)*(q+1); 
	
	
  
	double tol    = 100*DBL_EPSILON;
	
	if(fabs(xi[0]-knotU[numKnotU-1]) < tol) 
		xi[0] = knotU[numKnotU-1] - tol;
    
	if(fabs(xi[1]-knotV[numKnotV-1]) < tol) 
		xi[1] = knotV[numKnotV-1] - tol; 
	
	/* and evaluate the non-zero univariate B-spline basis functions
     * and first derivatives 
     */
	
	double *N      = (double *)malloc(sizeof(double)*(p+1));
    double *M      = (double *)malloc(sizeof(double)*(q+1));
	double **dersN = init2DArray(nU+1, p+1);
    double **dersM = init2DArray(nV+1, q+1);
	
	int spanU = FindSpan(nU, p, xi[0], knotU); 
    int spanV = FindSpan(nV, q, xi[1], knotV); 
    
	BasisFuns     (spanU, xi[0], p, knotU, N);
    BasisFuns     (spanV, xi[1], q, knotV, M);
    
	dersBasisFuns (spanU, xi[0], p, nU, knotU, dersN);	
    dersBasisFuns (spanV, xi[1], q, nV, knotV, dersM);
    
	
	/* and create NURBS approximation */
	
    int i, j, k, c;
    
    /*
    for(i=0;i<=p;i++){
        printf("dNdxi= %f\n", dersN[1][i]);
        printf("dNdet= %f\n", dersM[1][i]);
    }*/
    
    /*
    int uind = spanU - p;
    int vind;
    
    double w     = 0.0;
    double dwdxi = 0.0;
    double dwdet = 0.0;
    double wgt;
	
    for(j = 0; j <= q; j++)
    {
        vind = spanV - q + j;  
        
        for(i = 0; i <= p; i++)
        {               
            c   = uind + i + vind * (nU+1);            
            wgt = weight[c];
            
            w     += N[i]        * M[j] * wgt;
            dwdxi += dersN[1][i] * M[j] * wgt;
            dwdet += dersM[1][j] * N[i] * wgt;
        }
    }
    */
    
	/* create output */
    
    plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    
    double *R      = mxGetPr(plhs[0]);
    double *dRdxi  = mxGetPr(plhs[1]);
    double *dRdet  = mxGetPr(plhs[2]);
    
    k    = 0;
    
    double fac;
    
    /*printf("uind= %d\n", uind);  
    printf("vind= %d\n", vind);  
    printf("nU+1= %d\n", nU+1);  */
    
    for(j = 0; j <= q; j++)
    {
        for(i = 0; i <= p; i++)
        {               
            R[k]     = N[i]*M[j];
            dRdxi[k] = dersN[1][i]*M[j];
            dRdet[k] = dersM[1][j]*N[i];
            
            k += 1;
        }
    }
    
      	/*mexPrintf("\nWe have a knot vector with %d components\n", numWeights);
				for(k = 0; k < numWeights; k++) mexPrintf("%2.2f\t", weight[k]);
				mexPrintf("\n");*/
                
	free(N);
    free(M);
	free2Darray(dersN, (nU+1));
    free2Darray(dersM, (nV+1));	
}


 
 
 
 

