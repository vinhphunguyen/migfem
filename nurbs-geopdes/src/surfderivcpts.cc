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


DEFUN_DLD(surfderivcpts, args, nargout,"\
\nSURFDERIVCPTS: Compute control points of n-th derivatives of a NURBS surface.\n \
\n \
\nusage: pkl = surfderivcpts (n, p, U, m, q, V, P, d)  \
\n \
\n  INPUT:  \
\n\
\n        n+1, m+1 = number of control points \
\n        p, q     = spline order \
\n        U, V     = knots \
\n        P        = control points \
\n        d        = derivative order \
\n\
\n  OUTPUT: \
\n\
\n        pkl (k+1, l+1, i+1, j+1) = i,jth control point \
\n                                   of the surface differentiated k \
\n                                   times in the u direction and l \
\n                                   times in the v direction \
\n \
\n Adaptation of algorithm A3.7 from the NURBS book\n")
{
  //function pkl = surfderivcpts (n, p, U, m, q, V, P, d, r1, r2, s1, s2) 
  
  octave_value_list retval;

  octave_idx_type n = args(0).idx_type_value ();
  octave_idx_type p = args(1).idx_type_value ();
  RowVector U = args(2).row_vector_value (false, true);
  octave_idx_type m = args(3).idx_type_value ();
  octave_idx_type q = args(4).idx_type_value ();
  RowVector V = args(5).row_vector_value (false, true);
  Matrix P = args(6).matrix_value ();
  octave_idx_type d = args(7).idx_type_value ();
  
  octave_idx_type r1(0), r2 (n), s1 (0), s2 (m);
  
  
  if (args.length () == 12)
    {
      r1 = args (8).idx_type_value ();
      r2 = args (9).idx_type_value ();
      s1 = args (10).idx_type_value ();
      s2 = args (11).idx_type_value ();      
    } 
  else if  (args.length () > 8) 
    print_usage ();

  if (! error_state)
    {
      
      NDArray pkl;

      surfderivcpts (n, p, U, m, q, V, P, d, r1, r2, s1, s2,  pkl);

      retval(0) = octave_value (pkl);
    }
  return retval;
}


/*
%!test
%! plane = nrbdegelev(nrb4surf([0 0], [0 1], [1 0], [1 1]), [1, 1]);
%! 
%! pkl = surfderivcpts (plane.number(1)-1, plane.order(1)-1,
%!                       plane.knots{1}, plane.number(2)-1,
%!                     plane.order(2)-1, plane.knots{2}, 
%!                       squeeze (plane.coefs(1,:,:)), 2);
%! 
%! 
%! pkl2 = [  0   0   0   1   0   0   0   0   0   0   0   0   1   0 ...
%! 	0   0   0   0   0   0   0   1   0   0   0   0   0 0.5   0 ...
%! 	0   1   0   0   0   0   0 0.5   0   0   1   0   0   0   0 ...
%! 	0 0.5  0   0   1   0   0   0   0   0   1   0   0   0   0 ...
%! 	0   0   0   0   1   0   0   0   0   0   0   0   0   1   0 ...
%! 	0   0   0   0   0   0   0]';
%! 
%! assert (pkl(:),pkl2);
*/
