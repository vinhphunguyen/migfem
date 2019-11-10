#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "NURBS.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*	Return the NURBS basis function to matlab
	
	
	 We expect the function to be called as 
	 [interp interp_deriv] = NURBSinterpolation(xi, p, knot, points, weights)
	
	xi   = point where we want to interpolate
		knot = knot vector
		vector of points in format [pt1 pt2 pt3 pt4 .... ptn] 
		Robert Simpson, Cardiff University, UK */
	
	if(nrhs != 5) mexErrMsgTxt("You fool! You haven't passed in 3 arguments to the function."
		"We expect it to be in the form [interp interp_deriv] = NURBSinterpolation(xi, knot, points)\n");
			
	/* First get the inputs */
	
	double *xi   = mxGetPr(prhs[0]);	

	double *p_in = (mxGetPr(prhs[1]));
	int    p     = (int) *p_in;
	
	double *knot   = mxGetPr(prhs[2]);
	int    numKnot = mxGetN(prhs[2]);
	int    m       = numKnot - 1;
	int    n       = m - p -1;
	
	double *points   = mxGetPr(prhs[3]);
	int    numPoints = mxGetN(prhs[3]);
	
	double *weight   = mxGetPr(prhs[4]);
	int numWeights   = mxGetN(prhs[4]);
	
	double tol    = 100*DBL_EPSILON;
	
	if(fabs(*xi-knot[numKnot-1]) < tol) 
		*xi = knot[numKnot-1] - tol;
	
	/* and evaluate the non-zero basis functions*/
	
	double *N           = (double *)malloc(sizeof(double)*(p+1));
	double *NURBS       = (double *)malloc(sizeof(double)*(p+1));
	double *NURBS_deriv = (double *)malloc(sizeof(double)*(p+1));
	double **ders = init2DArray(n+1, p+1);
	
	int span = FindSpan(n, p, *xi, knot); 
	BasisFuns     (span, *xi, p, knot, N);
	dersBasisFuns (span, *xi, p, n, knot, ders);	
	
	/* and create NURBS approximation */
	int k, c;
	
	for(k = 0; k <=p; k++)
	{
		double w_interp = 0.0, dw_interp_dxi= 0.0;
		
		for(c = 0; c <= p; c++)
		{
			w_interp      += N[c] * weight[span-p+c];
			dw_interp_dxi += ders[1][c] * weight[span-p+c];
		}

		NURBS[k]       = N[k] * weight[span-p+k] / w_interp;
		NURBS_deriv[k] = weight[span-p+k] * ( w_interp * ders[1][k] - dw_interp_dxi * N[k] ) / (w_interp * w_interp);
	}

	double interp = 0.0;
	double interp_deriv = 0.0;
	
	for(c = 0; c <= p; c++)
	{
		interp       += NURBS[c] * points[span-p+c];
		interp_deriv += NURBS_deriv[c] * points[span-p+c];
	}
	
	free(N);
	free(NURBS);
	free(NURBS_deriv);
	free2Darray(ders, (n+1));
	
	plhs[0] = mxCreateDoubleScalar(interp); 
	plhs[1] = mxCreateDoubleScalar(interp_deriv);
}


 
 
 
 

