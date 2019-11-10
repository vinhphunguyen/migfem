#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "NURBS.h"

/*
 * Robert Simpson, Cardiff University, UK
 */

int main(int arg, char* argv[])
{
	int p = 2;																	//	our basis order
	double U[] = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0,5.0, 5.0};		//	knot vector
	int m = sizeof(U)/sizeof(double) - 1;										//	number of knot points
	int n = m - p - 1;															//	number of basis functions
	
	double u; 
	int c, span;
	
	/*printf("Number of knot points: %d\n", m);
	
	//
	//	Test our knot span index code
	//
	
	for(c = 0; c < m+1; ++c)
	{
		span = FindSpan(n, p, U[c], U);
		printf("For Xi=%2.1f, the knot span index is i=%d\n\n",U[c], span);
	}*/
	
	//
	//	Now let's test our routines which calculate the basis functions and derivatives
	//
	
	int numNonZeroBsFns = p + 1;
	
	/*
	double *N = (double*)malloc(sizeof(double) * numNonZeroBsFns);
	
	double **ders = init2DArray(n+1, p+1);
	
	for(c = 0; c < numNonZeroBsFns; ++c) N[c] = 0.0;
	
	double point[] = {0.0, 1.0, 2,0, 3.0, 4.0, 5.0};	//	some sample points to evaluate at
	int k = sizeof(point) / sizeof(double);				//	number of points (for our loop)
	int j,d;
	
	for(j = 0; j < k; ++j)
	{
		span = FindSpan(n, p, point[j], U);				//	get the span
		BasisFuns( span, point[j], p, U, N);			//	get the non-zero basis functions
		dersBasisFuns(span, point[j], p, n, U, ders);	//	get the derivatives ders[k][j]

		//	print out basis functions
		
		printf("The non-zero basis functions in the range %2.2f and %2.2f and u=%2.2f are:\n",U[span], U[span+1], point[j]);
		for(c = 0; c < numNonZeroBsFns; ++c) printf("\n %2.2f\n", N[c]);
		
		
		// 	print out derivatives
		
		for(d = 0; d < (n+1); d++)
		{
			printf("For k=%d, derivatives are: \n", d);
			for(c = 0; c < numNonZeroBsFns; c++) printf("%2.2f\t", ders[d][c]);
			printf("\n");
		}
		printf("\n\n");

	}
	
	free(N);
	free2Darray(ders, n+1);
	*/
	
	double *dN = (double*)malloc(sizeof(double) * (n+1));
	double **ders = init2DArray(n+1, p+1);
	
	double xi, Nip;
	int i=2;
	
	printf("xi\tN\tdNdxi\n");
	for(xi=0.0; xi <= 5.0; xi+=0.1)
	{
		span = FindSpan(n, p, xi, U);
		Nip = OneBasisFun(p, m, U, i, xi);
		dersOneBasisFuns(p, m, U, i, xi, n, dN);

	}
	
	xi = 2.1;
	span = FindSpan(n, p, xi, U);
	dersBasisFuns(span, xi, p, n, U, ders);	//	get the derivatives ders[k][j]
	printf("\n\n%2.2f\t%2.20f\t%2.20f\t%2.20f\n",xi, ders[1][0], ders[1][1], ders[1][2]);
		
	free(dN); 
	free2Darray(ders, n+1);
	return 0;
}


 
 
 
 
 

