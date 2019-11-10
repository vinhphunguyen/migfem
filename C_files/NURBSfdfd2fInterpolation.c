#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "NURBS.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*	Return the NURBS interpolations: c(u), c'(u) and c''(u) to matlab
	
	//
	// We expect the function to be called as 
	// [cu cu_der cu2_der] = NURBSfdfd2fInterpolation(xi, p, knot, points, weights)
	//
	//	xi   = point where we want to interpolate
	//	knot = knot vector
	//	vector of points in format [pt1 pt2 pt3 pt4 .... ptn] */
	
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
	
	double *N            = (double *)malloc(sizeof(double)*(p+1));
	double *NURBS        = (double *)malloc(sizeof(double)*(p+1));
	double *NURBS_deriv  = (double *)malloc(sizeof(double)*(p+1));
    double *NURBS_deriv2 = (double *)malloc(sizeof(double)*(p+1));
	double **ders = init2DArray(n+1, p+1);
	
	int span = FindSpan(n, p, *xi, knot); 
	BasisFuns     (span, *xi, p, knot, N);
	dersBasisFuns (span, *xi, p, n, knot, ders);	
	
	/* and create NURBS approximation */
	int k, c;
    double wc,wk,pc;
	
	for(k = 0; k <=p; k++)
	{
		double w_interp       = 0.0;
        double dw_interp_dxi  = 0.0;
        double dw2_interp_dxi = 0.0;
		
		for(c = 0; c <= p; c++)
		{
            wc              = weight[span-p+c];
			w_interp       += N[c]       * wc;
			dw_interp_dxi  += ders[1][c] * wc;
            dw2_interp_dxi += ders[2][c] * wc;
		}

        wk              = weight[span-p+k];
		NURBS[k]        = N[k] * wk / w_interp;
		NURBS_deriv[k]  = wk * ( w_interp * ders[1][k] - dw_interp_dxi * N[k] ) / (w_interp * w_interp);
        NURBS_deriv2[k] = ( wk*ders[2][k] - 2*dw_interp_dxi * NURBS_deriv[k] - dw2_interp_dxi*NURBS[k]) / w_interp;
	}

	double interp        = 0.0;
	double interp_deriv  = 0.0;
    double interp_deriv2 = 0.0;
	
	for(c = 0; c <= p; c++)
	{
        pc             = points[span-p+c];
		interp        += NURBS[c]        * pc;
		interp_deriv  += NURBS_deriv[c]  * pc;
        interp_deriv2 += NURBS_deriv2[c] * pc;
	}
	
	free(N);
	free(NURBS);
	free(NURBS_deriv);
    free(NURBS_deriv2);
	free2Darray(ders, (n+1));
	
	plhs[0] = mxCreateDoubleScalar(interp); 
	plhs[1] = mxCreateDoubleScalar(interp_deriv);
    plhs[2] = mxCreateDoubleScalar(interp_deriv2);
}


 
 
 
 

