#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "NURBS.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*	Return the NURBS basis function to matlab
	
	//
	// We expect the function to be called as 
	// [interp] = NURBSinterpolation(xi, eta, p, q, uknot, vknot
                                   pointsX, pointsY, weights)
	//
	//	xi   = point where we want to interpolate [xi eta]
	//	uknot = knot vector u direction
    //	vknot = knot vector v direction
	//	vector of points in format [pt1 pt2 pt3 pt4 .... ptn] */
	
	if(nrhs != 9) mexErrMsgTxt("You fool! You haven't passed in 9 arguments to the function."
		"We expect it to be in the form [interp interp_deriv]=NURBSinterpolation(xi, p, q, uknot, vknot, points, weights)\n");
			
	/* First get the inputs */
	
	double *xi      = mxGetPr(prhs[0]);	
    double *eta     = mxGetPr(prhs[1]);	
	double *p_in    = mxGetPr(prhs[2]);
    double *q_in    = mxGetPr(prhs[3]);
    double *uknot   = mxGetPr(prhs[4]);
    double *vknot   = mxGetPr(prhs[5]);
    double *pointsX = mxGetPr(prhs[6]);
    double *pointsY = mxGetPr(prhs[7]);
    double *weight  = mxGetPr(prhs[8]);
    
    /* create output */
    
    plhs[0] = mxCreateDoubleMatrix(1,2,mxREAL); 
    
    double *interp  = mxGetPr(plhs[0]);
    
    /* Then use them as in standard C*/
    
	int    p       = (int) *p_in;
    int    q       = (int) *q_in;
    
	int    numKnotU = mxGetN(prhs[4]);
    int    numKnotV = mxGetN(prhs[5]);
    
	int    mu       = numKnotU - 1;
    int    mv       = numKnotV - 1;
    
	int    nu       = mu - p - 1;
    int    nv       = mv - q - 1;
	
	
	int    numPoints  = mxGetN(prhs[6]);
	int    numWeights = mxGetN(prhs[7]);
    
	
	double tol        = 100*DBL_EPSILON;
	
	if(fabs(*xi -uknot[numKnotU-1]) < tol) *xi  = uknot[numKnotU-1] - tol;
    if(fabs(*eta-vknot[numKnotV-1]) < tol) *eta = vknot[numKnotV-1] - tol;
	
	/* and evaluate the non-zero B-spline basis functions*/
	
	double *N        = (double *)malloc(sizeof(double)*(p+1));
    double *M        = (double *)malloc(sizeof(double)*(q+1));
	
	int spanU = FindSpan(nu, p, *xi,  uknot); 
    int spanV = FindSpan(nv, q, *eta, vknot); 
    
	BasisFuns     (spanU, *xi,  p, uknot, N);
    BasisFuns     (spanV, *eta, q, vknot, M);
    
    /* and compute the approximation */
    
	interp[0] = 0.0;
    interp[1] = 0.0;
    
    double wght      = 0.0;
    
    double tempu, tempv, tempw;
	
    int    vind, uind = spanU - p;;
    int    k,l,id;
    
	for(l = 0; l <= q; l++)
    {
        tempu = 0;
        tempv = 0;
        tempw = 0;
 
        vind  = spanV - q + l;
        
        //printf("vind = %d\n", vind);
        //printf("uind = %d\n", uind);
        
        for(k = 0; k <= p; k++)
	    {
           //access control point P(vind,uind+k)
		   
           id      = (uind+k)*(nv+1) + vind+1; 
           tempu  += N[k] * (*pointsX+vind*nv+uind+k) * weight[id];
           tempv  += N[k] * (*pointsY+vind*nv+uind+k) * weight[id];
           tempw  += N[k] * weight[id];
	    }
        
        interp[0]  += tempu * M[l];
        interp[1]  += tempv * M[l];
        wght       += tempw * M[l];
    }
    
    // projection
            
    interp[0] /= wght;
    interp[1] /= wght;
	
	free(N);
    free(M);
}


 
 
 
 

