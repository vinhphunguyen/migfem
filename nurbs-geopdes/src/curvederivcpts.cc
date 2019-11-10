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


DEFUN_DLD(curvederivcpts, args, nargout,"\
\nCURVEDERIVCPTS: Compute control points of n-th derivatives of a B-spline curve.\n \
\n \
\n usage: pk = curvederivcpts (n, p, U, P, d) \
\n        pk = curvederivcpts (n, p, U, P, d, r1 r2) \
\n \
\n If r1, r2 are not given, all the control points are computed. \
\n \
\n  INPUT: \
\n         n+1 = number of control points \
\n         p   = degree of the spline \
\n         d   = maximum derivative order (d<=p) \
\n         U   = knots \
\n         P   = control points \
\n         r1  = first control point to compute \
\n         r2  = auxiliary index for the last control point to compute \
\n\
\n  OUTPUT: \
\n         pk(k,i) = i-th control point (k-1)-th derivative, r1 <= i <= r2-k \
\n \
\n Adaptation of algorithm A3.3 from the NURBS book\n")

{
  
  octave_value_list retval;

  octave_idx_type n = args(0).idx_type_value ();
  octave_idx_type p = args(1).idx_type_value ();
  RowVector U = args(2).row_vector_value (false, true);
  NDArray P = args(3).array_value ();
  octave_idx_type d = args(4).idx_type_value ();

  octave_idx_type r1(0), r2(n);
  if (args.length () == 7)
    {
      r1 = args (5).idx_type_value ();
      r2 = args (6).idx_type_value ();
    }
  else  if (args.length () > 5)
    print_usage ();

  if (! error_state)  
    {
      octave_idx_type r = r2 - r1;
      Matrix pk (d+1 <= r+1 ? d+1 : r+1, r+1, 0.0);

      curvederivcpts (n, p, U, P, d, r1, r2, pk);

      retval(0) = octave_value (pk);
    }
  return retval;
}

/*
%!test
%! line = nrbmak([0.0 1.5; 0.0 3.0],[0.0 0.0 1.0 1.0]);
%! pk   = curvederivcpts (line.number-1, line.order-1, line.knots,
%!                        line.coefs(1,:), 2);
%! assert (pk, [0 3/2; 3/2 0], 100*eps);
*/
