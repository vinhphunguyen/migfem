/* Copyright (C) 2009 Carlo de Falco
   some functions are adapted from the m-file implementation which is
   Copyright (C) 2000 Mark Spink, 2007 Daniel Claxton

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <octave/oct.h>
#include "low_level_functions.h"
#include <iostream>

octave_idx_type findspan(int n, int p, double u, const RowVector& U)

// Find the knot span of the parametric point u. 
//
// INPUT:
//
//   n - number of control points - 1
//   p - spline degree       
//   u - parametric point    
//   U - knot sequence
//
// RETURN:
//
//   s - knot span
//
// Note: This is NOT
// Algorithm A2.1 from 'The NURBS BOOK' pg68
// as that algorithm only works for nonperiodic
// knot vectors, nonetheless the results should 
// be EXACTLY the same if U is nonperiodic

/*
Below is the original implementation from the NURBS Book
{
  int low, high, mid;
  // special case
  if (u == U(n+1)) return(n);

  // do binary search
  low = p;
  high = n + 1;
  mid = (low + high) / 2;
  while (u < U(mid) || u >= U(mid+1))
    {

      if (u < U(mid))
	high = mid;
      else
	low = mid;
      mid = (low + high) / 2;
    }  

  return(mid);
}
*/

{
  // FIXME : this implementation has linear, rather than log complexity
  int ret = 0;
  while ((ret++ < n) && (U(ret) <= u)) {
  };
  return (ret-1);
}

void basisfun(int i, double u, int p, const RowVector& U, RowVector& N)

// Basis Function. 
//
// INPUT:
//
//   i - knot span  ( from FindSpan() )
//   u - parametric point
//   p - spline degree
//   U - knot sequence
//
// OUTPUT:
//
//   N - Basis functions vector[p+1]
//
// Algorithm A2.2 from 'The NURBS BOOK' pg70.
{
  int j,r;
  double saved, temp;

  // work space
  OCTAVE_LOCAL_BUFFER(double, left,  p+1);
  OCTAVE_LOCAL_BUFFER(double, right, p+1);
  
  N(0) = 1.0;
  for (j = 1; j <= p; j++)
    {
      left[j]  = u - U(i+1-j);
      right[j] = U(i+j) - u;
      saved = 0.0;
      
      for (r = 0; r < j; r++)
	{
	  temp = N(r) / (right[r+1] + left[j-r]);
	  N(r) = saved + right[r+1] * temp;
	  saved = left[j-r] * temp;
	} 
      
      N(j) = saved;
    }

}


void basisfunder (int i, int pl, double u, const RowVector& u_knotl, int nders, NDArray& ders)
{
  
//    BASISFUNDER:  B-Spline Basis function derivatives
//
//     INPUT:
//
//       i   - knot span
//       pl  - degree of curve
//       u   - parametric points
//       k   - knot vector
//       nd  - number of derivatives to compute
//
//     OUTPUT:
//
//       ders - ders(n, i, :) (i-1)-th derivative at n-th point


  //     ders = zeros(nders+1,pl+1);
  Matrix ndu(octave_idx_type(pl+1), octave_idx_type(pl+1), 0.0); // ndu = zeros(pl+1,pl+1);
  RowVector left(octave_idx_type(pl+1), 0.0);                    // left = zeros(pl+1);
  RowVector right(left);                                         // right = zeros(pl+1);
  Matrix a(2, octave_idx_type(pl+1), 0.0);                       // a = zeros(2,pl+1);
  double saved = 0.0, d = 0.0, temp = 1.0;
  octave_idx_type s1(0), s2(1), rk, pk, j, k, r, j1, j2;

  
  ndu(0,0) = 1;                                                   // ndu(1,1) = 1;
  
  for (j=1; j<=pl; j++)                                           // for j = 1:pl
    {
      left(j) = u - u_knotl(i+1-j);                               // left(j+1) = u - u_knotl(i+1-j);
      right(j) = u_knotl(i+j) - u;                                // right(j+1) = u_knotl(i+j) - u;
      saved = 0.0;                                                // saved = 0;
      for (r=0; r<=j-1; r++)                                      // for r = 0:j-1
	{
	  ndu(j, r) = right(r+1) + left(j-r);                     // ndu(j+1,r+1) = right(r+2) + left(j-r+1);
	  temp = ndu(r,j-1)/ndu(j,r);                             // temp = ndu(r+1,j)/ndu(j+1,r+1);
	  ndu(r,j) = saved + right(r+1)*temp;                     // ndu(r+1,j+1) = saved + right(r+2)*temp;
	  saved = left(j-r)*temp;                                 // saved = left(j-r+1)*temp;
	}                                                         // end
      ndu(j,j) = saved;                                           // ndu(j+1,j+1) = saved;
    }                                                             // end

  for (j=0; j<=pl; j++)                                           // for j = 0:pl
    ders(0,j) = ndu(j,pl);                                        // ders(1,j+1) = ndu(j+1,pl+1);
                                                                  // end
      
  for (r=0; r<=pl; r++)                                           // for r = 0:pl
    {
      s1 = 0;                                                     // s1 = 0;
      s2 = 1;                                                     // s2 = 1;
      a(0,0) = 1;                                                 // a(1,1) = 1;
        for (k=1; k<=nders; k++)                                  // for k = 1:nders %compute kth derivative
	  {

	    d = 0.0;                                              // d = 0;
	    rk = r-k;                                             // rk = r-k;
	    pk = pl - k;                                          // pk = pl-k;
	    
	    if (r >= k)                                           // if (r >= k)
	      {
		a(s2, 0) = a(s1, 0)/ndu(pk+1,rk);                 // a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1);
		d = a(s2, 0)*ndu(rk,pk);                          // d = a(s2+1,1)*ndu(rk+1,pk+1);
	      }                                                   // end
	    
	    if (rk >= -1)                                         // if (rk >= -1)
	      j1 = 1;                                             // j1 = 1;
	    else                                                  // else 
	      j1 = -rk;                                           // j1 = -rk;
	                                                          // end

	    if ((r-1) <= pk)                                      // if ((r-1) <= pk)
	      j2 = k-1;                                           // j2 = k-1;
	    else                                                  // else 
	     j2 = pl-r;                                           // j2 = pl-r;
                                                                  // end

	    for (j=j1; j <= j2; j++)                              // for j = j1:j2
	      {
		a(s2,j) = (a(s1,j) - a(s1,j-1))/ndu(pk+1,rk+j);   // a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j))/ndu(pk+2,rk+j+1);
		d += a(s2,j)*ndu(rk+j,pk);                        // d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1);
	      }                                                   // end

	    if (r <= pk)                                          // if (r <= pk)
	      {
		a(s2,k) = -a(s1,k-1)/ndu(pk+1,r);                 // a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1);
		d += a(s2,k)*ndu(r,pk);                           // d = d + a(s2+1,k+1)*ndu(r+1,pk+1);
	      }	                                                  // end
	    
	    ders(k,r) = d;                                        // ders(k+1,r+1) = d;
	    j = s1;                                               // j = s1;
	    s1 = s2;                                              // s1 = s2;
	    s2 = j;                                               // s2 = j;
	  }                                                       // end
    }                                                             // end

  r = pl;                                                         // r = pl;
  for (k=1; k <= nders; k++)                                      // for k = 1:nders
    {
      for (j=0; j<=pl; j++)                                       // for j = 0:pl
	ders(k,j) = ders(k,j)*r;                                  // ders(k+1,j+1) = ders(k+1,j+1)*r;
                                                                  // end
      r = r*(pl-k);                                               // r = r*(pl-k);
    }                                                             // end
  
}


int curvederivcpts (octave_idx_type n, octave_idx_type p, 
		    const RowVector &U, const NDArray &P, 
		    octave_idx_type  d, octave_idx_type r1, 
		    octave_idx_type r2, Matrix &pk)
{
  octave_idx_type r = r2 - r1;
  for (octave_idx_type i(0); i<=r; i++)
    pk(0, i) = P(r1+i);

  for (octave_idx_type k (1); k<=d; k++)
    {
      octave_idx_type tmp = p - k + 1;
      for (octave_idx_type i (0); i<=r-k; i++)
	{
	  pk (k, i) = tmp * (pk(k-1,i+1)-pk(k-1,i)) / 
	    (U(r1+i+p+1)-U(r1+i+k));
	}
    }
  return 0;
}



int surfderivcpts (octave_idx_type n, octave_idx_type  p, const RowVector& U, 
		   octave_idx_type m, octave_idx_type q, const RowVector& V, 
		   const Matrix& P, octave_idx_type d, octave_idx_type r1, 
		   octave_idx_type r2, octave_idx_type s1, 
		   octave_idx_type s2, NDArray &pkl)
{
  octave_idx_type r = r2-r1, s = s2-s1;  
  
  octave_idx_type du = d <= p ? d : p;   
  octave_idx_type dv = d <= q ? d : q; 
  
  dim_vector idxa (4, 1);

  Array<octave_idx_type> idxta (idxa, 0);
  Array<idx_vector> idxva (idxa, idx_vector (':'));
  
  idxa.resize (4);
  idxa(0) = (du+1); idxa(1) = (dv+1); 
  idxa(2) = (r+1);  idxa(3) = (s+1); 
  
  pkl.resize (idxa, 0.0); 
      
  for (octave_idx_type j(s1); j<=s2; j++)
    {
    
      Matrix temp (du <= n ? (du+1) : (n+1), n+1, 0.0);
      curvederivcpts (n, p, U, P.extract (0, j, P.rows()-1, P.cols ()-1), du, r1, r2, temp);
      for (octave_idx_type k(0); k<=du; k++)
	{

	  for ( octave_idx_type i(0); i<=r-k; i++)
	    {
	      assert (k<idxa(0) && i<idxa(2) && j-s1<idxa(3));
	      idxta (0) = k; idxta (1) = 0;
	      idxta (2) = i; idxta (3) = j-s1;
	      pkl(idxta) = temp (k, i);
	    }
	}
    }

  for (octave_idx_type k (0); k<=du; k++)
    {

      for (octave_idx_type i(0); i<=r-k; i++)
	{

	  octave_idx_type dd = (d-k) <= dv ? (d-k) : dv;
	  Matrix temp (dd <= m ? (dd+1) : (m+1), m+1, 0.0);
	  
	  idxva (0) = idx_vector(k); idxva (1) = idx_vector(0);
	  idxva (2) = idx_vector(i); idxva (3) = idx_vector(':');
	  NDArray temp2 (pkl.index (idxva));
	  curvederivcpts (m, q, V.extract (s1, V.length ()-1), temp2.squeeze (), dd, 0, s, temp);
	  
	  for (octave_idx_type l(1); l<=dd; l++)
	    {

	      for (octave_idx_type j(0); j<=s-l; j++)
		{
		  assert (k<idxa(0) && l<idxa(1) && i<idxa(2) && j<idxa(3));
		  idxta (0) = k; idxta (1) = l;
		  idxta (2) = i; idxta (3) = j;
		  pkl(idxta) = temp (l, j);
		}
	    }
	}
    }
  return (0);
}


int surfderiveval (octave_idx_type n, octave_idx_type p, const RowVector &U, 
		   octave_idx_type m, octave_idx_type q, const RowVector &V, 
		   const Matrix &P, double u, double v, octave_idx_type d, 
		   Matrix &skl)
{

  
  Array<octave_idx_type> idx(dim_vector (4, 1), 0);
  octave_idx_type du = d <= p ? d: p;   
  octave_idx_type dv = d <= q ? d: q;
  
  skl.resize (d+1, d+1, 0.0);
  
  octave_idx_type uspan = findspan (n, p, u, U);
  Matrix Nu (p+1, p+1, 0.0);
  for (octave_idx_type ip(0); ip<=p; ip++)
    {
      RowVector temp (ip+1, 0.0);
      basisfun (uspan, u, ip, U, temp);
      Nu.insert (temp.transpose (), 0,ip);
    }
  
  octave_idx_type vspan = findspan (m, q, v, V);
  Matrix Nv (q+1, q+1, 0.0);
  for (octave_idx_type iq(0); iq<=q; iq++)
    {
      RowVector temp (iq+1, 0.0);
      basisfun (vspan, v, iq, V, temp);
      Nv.insert (temp.transpose (), 0, iq);
    }
  NDArray pkl;
  surfderivcpts (n, p, U, m, q, V, P, d, uspan-p, uspan, vspan-q, vspan, pkl);
  
  for (octave_idx_type k(0); k<=du; k++)
    {
      octave_idx_type dd = d-k <= dv ? d-k : dv;
      for (octave_idx_type l(0);l <= dd; l++)
	{
	  skl(k,l) = 0.0;
	  for (octave_idx_type i(0); i<=q-l; i++)
	    {
	      double tmp = 0.0;
	      for (octave_idx_type j(0); j<=p-k; j++)
		{
		  idx(0) = k; idx(1)=l; idx(2)=j; idx(3) =i;
		  tmp += Nu(j,p-k) * pkl(idx);
		}
	      skl(k,l) += Nv(i,q-l) * tmp;
	    }
	}
    }
  return (0);
}
