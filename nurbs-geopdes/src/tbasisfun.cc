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
#include <iostream>

double onebasisfun__ (double u, octave_idx_type p, RowVector U)
{

  //std::cout << "u=" << u << " " << "p=" << p << " \n" << "U=" << U;
  
  double N = 0.0;
  if ((u < U.min ()) || ( u > U.max ()))
    return (N);
  else if (p == 0)
    return (1.0);
 
  double ln = u - U(0);
  double ld = U(U.length () - 2) - U(0);
  if (ld != 0)
    N += ln * onebasisfun__ (u, p-1, U.extract (0, U.length () - 2))/ ld; 
    
  double dn = U(U.length () - 1) - u;
  double dd = U(U.length () - 1) - U(1);
  if (dd != 0)
    N += dn * onebasisfun__ (u, p-1, U.extract (1, U.length () - 1))/ dd;
    
  return (N);
}

   
DEFUN_DLD(tbasisfun, args, nargout,"\
TBASISFUN: Compute a B- or T-Spline basis function from its local knot vector.\n\
\n\
 usage:\n\
\n\
 N = tbasisfun (u, p, U)\n\
 N = tbasisfun ([u; v], [p q], {U, V})\n\
 \n\
 INPUT:\n\
  u or [u; v] : points in parameter space where the basis function is to be\n\
  evaluated \n\
  \n\
  U or {U, V} : local knot vector\n\
\n\
 p or [p q] : polynomial order of the basis function\n\
\n\
 OUTPUT:\n\
  N : basis function evaluated at the given parametric points\n")

{
  
  octave_value_list retval;
  Matrix u = args(0).matrix_value ();

  if (! args(2).is_cell ())
    {

      double p = args(1).idx_type_value ();
      RowVector U = args(2).row_vector_value (true, true);
      assert (U.numel () == p+2);
      
      RowVector N(u.cols ());
      for (octave_idx_type ii=0; ii<u.numel (); ii++)
	N(ii) = onebasisfun__ (u(ii), p, U);

    }  else {

    RowVector p = args(1).row_vector_value ();
    Cell C = args(2).cell_value ();
    RowVector U = C(0).row_vector_value (true, true);
    RowVector V = C(1).row_vector_value (true, true);
    

    RowVector N(u.cols ());
    for (octave_idx_type ii=0; ii<u.cols (); ii++)
      {
	N(ii) = onebasisfun__ (u(0, ii), octave_idx_type(p(0)), U) *
	  onebasisfun__ (u(1, ii), octave_idx_type(p(1)), V);
	//std::cout << "N=" << N(ii) << "\n\n\n";
      }
    retval(0) = octave_value (N);
  }
  return retval;
}


/*
%!demo
%! U = {[0 0 1/2 1 1], [0 0 0 1 1]};
%! p = [3, 3];
%! [X, Y] = meshgrid (linspace(0, 1, 30));
%! u = [X(:), Y(:)]';
%! N = tbasisfun (u, p, U);
%! surf (X, Y, reshape (N, size(X)))
*/
