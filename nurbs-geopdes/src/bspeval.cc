/* Copyright (C) 2009 Carlo de Falco
  
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
#include <omp.h>

static bool bspeval_bad_arguments(const octave_value_list& args);

DEFUN_DLD(bspeval, args, nargout,"\
 BSPEVAL:  Evaluate B-Spline at parametric points\n\
\n\
\n\
 Calling Sequence:\n\
\n\
   p = bspeval(d,c,k,u)\n\
\n\
    INPUT:\n\
\n\
       d - Degree of the B-Spline.\n\
       c - Control Points, matrix of size (dim,nc).\n\
       k - Knot sequence, row vector of size nk.\n\
       u - Parametric evaluation points, row vector of size nu.\n\
 \n\
    OUTPUT:\n\
\n\
       p - Evaluated points, matrix of size (dim,nu)\n\
")
{

  octave_value_list retval;
  if (!bspeval_bad_arguments (args))
    {      
      int             d = args(0).int_value();
      const Matrix    c = args(1).matrix_value();
      const RowVector k = args(2).row_vector_value();
      const NDArray   u = args(3).array_value();
      
      octave_idx_type nu = u.length();
      octave_idx_type mc = c.rows(),
        nc = c.cols();
      
      Matrix p(mc, nu, 0.0);
      
      if (!error_state)
        {
          if (nc + d == k.length() - 1) 
            {	 
#pragma omp parallel default (none) shared (d, c, k, u, nu, mc, nc, p)
              {
                RowVector N(d+1,0.0);
                int s, tmp1;
                double tmp2;
#pragma omp for 
                for (octave_idx_type col=0; col<nu; col++)
                  {	
                    //printf ("thread %d, col %d\n", omp_get_thread_num (), col);
                    s = findspan (nc-1, d, u(col), k);
                    basisfun (s, u(col), d, k, N);    
                    tmp1 = s - d;                
                    for (octave_idx_type row(0); row<mc; row++)
                      {
                        tmp2 = 0.0;
                        for ( octave_idx_type i(0); i<=d; i++)                   
                          tmp2 +=  N(i)*c(row,tmp1+i);	  
                        p(row,col) = tmp2;
                      }             
                  }   
              }// end omp 
            }
          else 
            {
              error("inconsistent bspline data, d + columns(c) != length(k) - 1.");
            }
          retval(0) = octave_value(p);
        }
    }      
  return retval;
} 

static bool bspeval_bad_arguments (const octave_value_list& args) 
{ 
  if (args.length() != 4)
    {
      error("bspeval: wrong number of input arguments.");
      return true;
    }
  if (!args(0).is_real_scalar()) 
    { 
      error("bspeval: degree should be a scalar."); 
      return true; 
    } 
  if (!args(1).is_real_matrix()) 
    { 
      error("bspeval: the control net should be a matrix of doubles."); 
      return true; 
    } 
  if (!args(2).is_real_matrix()) 
    { 
      error("bspeval: the knot vector should be a real vector."); 
      return true; 
    } 
  if (!args(3).is_real_type()) 
    { 
      error("bspeval: the set of parametric points should be an array of doubles."); 
      return true; 
    } 
  return false; 
} 
