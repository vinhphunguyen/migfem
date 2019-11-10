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
#include <iostream>


DEFUN_DLD(surfderiveval, args, nargout,"\
\nSURFDERIVEVAL: Compute the derivatives of a B-spline surface\
\n\
\n usage: skl = surfderiveval (n, p, U, m, q, V, P, u, v, d) \
\n\
\n  INPUT: \
\n\
\n        n+1, m+1 = number of control points\
\n        p, q     = spline order\
\n        U, V     = knots\
\n        P        = control points\
\n        u,v      = evaluation points\
\n        d        = derivative order\
\n\
\n  OUTPUT:\
\n\
\n        skl (k+1, l+1) =  surface differentiated k\
\n                          times in the u direction and l\
\n                          times in the v direction\
\n\
\n Adaptation of algorithm A3.8 from the NURBS book\n")
{

  //function skl = surfderiveval (n, p, U, m, q, V, P, u, v, d) 
  
  octave_value_list retval;

  octave_idx_type n = args(0).idx_type_value ();
  octave_idx_type p = args(1).idx_type_value ();
  RowVector U = args(2).row_vector_value (false, true);
  octave_idx_type m = args(3).idx_type_value ();
  octave_idx_type q = args(4).idx_type_value ();
  RowVector V = args(5).row_vector_value (false, true);
  Matrix P = args(6).matrix_value ();
  double u = args(7).double_value ();
  double v = args(8).double_value ();
  octave_idx_type d = args(9).idx_type_value ();

  if (! error_state)
    {
      Matrix skl;
      surfderiveval (n, p, U, m, q, V, P, u, v, d, skl);
      retval(0) = octave_value (skl);
    }
  return retval;
}

/*
%!shared srf
%!test
%! k = [0 0 0 1 1 1];
%! c = [0 1/2 1];
%! [coef(2,:,:), coef(1,:,:)] = meshgrid (c, c);
%! srf = nrbmak (coef, {k, k});
%! skl = surfderiveval (srf.number(1)-1, 
%!                      srf.order(1)-1, 
%!                      srf.knots{1}, 
%!                      srf.number(2)-1, 
%!                      srf.order(2)-1, 
%!                      srf.knots{2},
%!                      squeeze(srf.coefs(1,:,:)), .5, .5, 1) ;
%! assert (skl, [.5 0; 1 0])
%!test
%! srf = nrbkntins (srf, {[], rand(1,2)});
%! skl = surfderiveval (srf.number(1)-1, 
%!                      srf.order(1)-1, 
%!                      srf.knots{1},
%!                      srf.number(2)-1, 
%!                      srf.order(2)-1, 
%!                      srf.knots{2},
%!                      squeeze(srf.coefs(1,:,:)), .5, .5, 1) ;
%! assert (skl, [.5 0; 1 0], 100*eps)
*/
