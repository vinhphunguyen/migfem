#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "NURBS.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Return the 2D NURBS basis functions and first/second derivatives to matlab
	   All non-zero basis functions and derivatives at point [xi,eta] are computed.
	
	  We expect the function to be called as 
	  [R dRdxi dRdeta dR2dxi dR2deta dR2dxideta] = NURBS2DBasis2ndDers(...
                                xi, p, q, knotU,knotV, weights)
	
		xi           = point, [xi eta], where we want to interpolate
		knotU, knotV = knot vectors
		weights      = vector of weights 
    
     Vinh Phu Nguyen, nvinhphu@gmail.com
    */
	
	if(nrhs != 6) mexErrMsgTxt("You fool! You haven't passed in 6 arguments to the function."
		"We expect it to be in the form [dRdxi dRdeta] = NURBSinterpolation(xi, p, q, knotU, knotV, weights)\n");
			
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
	
	double *weight  = mxGetPr(prhs[5]);
    int    numWeights = mxGetN(prhs[5]);
		
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
	
	int spanU      = FindSpan(nU, p, xi[0], knotU); 
    int spanV      = FindSpan(nV, q, xi[1], knotV); 
    
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
    
    int uind = spanU - p;
    int vind;
    
    double w      = 0.0; /* w = N_I w_I*/
    double dwdxi  = 0.0; /* first derivative of w w.r.t xi*/
    double d2wdxi = 0.0; /* second derivative of w w.r.t xi*/
    double dwdet  = 0.0; /* first derivative of w w.r.t eta*/
    double d2wdet = 0.0; /* second derivative of w w.r.t eta*/
    double d2wdxe = 0.0; /* second derivative of w w.r.t xi-eta*/
    double wi;
	
    for(j = 0; j <= q; j++)
    {
        vind = spanV - q + j;  
        
        for(i = 0; i <= p; i++)
        {               
            c   = uind + i + vind * (nU+1);            
            wi  = weight[c];
            
            w      += N[i]        * M[j] * wi;
            dwdxi  += dersN[1][i] * M[j] * wi;
            d2wdxi += dersN[2][i] * M[j] * wi;
            dwdet  += dersM[1][j] * N[i] * wi;
            d2wdet += dersM[2][j] * N[i] * wi;
            d2wdxe += dersN[1][i] * dersM[1][j] * wi;
        }
    }

	/* create output */
    
    plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    plhs[3] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    plhs[4] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    plhs[5] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    
    double *R      = mxGetPr(plhs[0]);
    double *dRdxi  = mxGetPr(plhs[1]); /* first derivative to xi*/
    double *dRdet  = mxGetPr(plhs[2]); /* first derivative to eta*/
    double *dR2dxi = mxGetPr(plhs[3]);
    double *dR2det = mxGetPr(plhs[4]);
    double *dR2dxe = mxGetPr(plhs[5]);
    
    uind = spanU - p;
    k    = 0;
    
    double fac;
    
    /*printf("uind= %d\n", uind);  
    printf("vind= %d\n", vind);  
    printf("nU+1= %d\n", nU+1);  */
    
    double NiMj, w3Inv, w2Inv;
    
    for(j = 0; j <= q; j++)
    {
        vind = spanV - q + j; 
        
        for(i = 0; i <= p; i++)
        {               
            c        = uind + i + vind*(nU+1);            
            wi       = weight[c];
            
            NiMj     = N[i]*M[j];
            w3Inv    = 1/w/w/w;
            w2Inv    = w3Inv*w;
            fac      = weight[c]*w2Inv;
            
            R[k]     = NiMj*fac*w;
            
            dRdxi[k] = (dersN[1][i]*M[j]*w - NiMj*dwdxi) * fac;
            dRdet[k] = (dersM[1][j]*N[i]*w - NiMj*dwdet) * fac;
            
            dR2dxi[k] = wi*(dersN[2][i]*M[j]/w - 2*dersN[1][i]*M[j]*dwdxi*w2Inv - NiMj*d2wdxi*w2Inv + 2*NiMj*dwdxi*dwdxi*w3Inv);
            dR2det[k] = wi*(dersM[2][j]*N[i]/w - 2*dersM[1][j]*N[i]*dwdet*w2Inv - NiMj*d2wdet*w2Inv + 2*NiMj*dwdet*dwdet*w3Inv);
            dR2dxe[k] = wi*(dersN[1][i]*dersM[1][j]/w - dersN[1][i]*M[j]*dwdet*w2Inv - N[i]*dersM[1][j]*dwdxi*w2Inv - NiMj*d2wdxe*w2Inv + 2*NiMj*dwdxi*dwdet*w3Inv);
            
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


 
 
 
 

