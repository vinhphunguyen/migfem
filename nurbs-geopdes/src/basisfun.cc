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

DEFUN_DLD(basisfun, args, nargout, "\n\
 BASISFUN: Compute B-Spline Basis Functions \n\
\n\
Calling Sequence:\n\
\n\
  N = basisfun(iv,uv,p,U)\n\
\n\
 INPUT:\n\
\n\
   iv - knot span  ( from FindSpan() )\n\
   uv - parametric point\n\
   p - spline degree\n\
   U - knot sequence\n\
\n\
 OUTPUT:\n\
\n\
   N - Basis functions vector(numel(uv)*(p+1))\n\
\n\
 Algorithm A2.2 from 'The NURBS BOOK' pg70.\n\
\n\
")
{

  octave_value_list retval;
  const NDArray   i = args(0).array_value();
  const NDArray   u = args(1).array_value();
  int       p = args(2).idx_type_value();
  const RowVector U = args(3).row_vector_value();
  RowVector N(p+1, 0.0);
  Matrix    B(u.length(), p+1, 0.0);
  
  if (!error_state)
    {
      for (octave_idx_type ii(0); ii < u.length(); ii++)
	{
	  basisfun(int(i(ii)), u(ii), p, U, N);
	  B.insert(N, ii, 0);
	}
      
      retval(0) = octave_value(B);
    }
  return retval;
} 

/*
%!shared n, U, p, u, s
%!test
%!  n = 3; 
%!  U = [0 0 0 1/2 1 1 1]; 
%!  p = 2; 
%!  u = linspace(0, 1, 10);  
%!  s = findspan(n, p, u, U); 
%!  assert (s, [2*ones(1, 5) 3*ones(1, 5)]);
%!test
%!  Bref = [1.00000   0.00000   0.00000
%!          0.60494   0.37037   0.02469
%!          0.30864   0.59259   0.09877
%!          0.11111   0.66667   0.22222
%!          0.01235   0.59259   0.39506
%!          0.39506   0.59259   0.01235
%!          0.22222   0.66667   0.11111
%!          0.09877   0.59259   0.30864
%!          0.02469   0.37037   0.60494
%!          0.00000   0.00000   1.00000];
%!  B = basisfun(s, u, p, U);
%!  assert (B, Bref, 1e-5);
*/
