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
#include <octave/oct-map.h>
#include <octave/parse.h>
#include "low_level_functions.h"

DEFUN_DLD(nrb_srf_basisfun__, args, nargout,"\
 NRB_SRF_BASISFUN__:  Undocumented private function\
")
{

  octave_value_list retval, newargs;

  const NDArray points = args(0).array_value();
  const Octave_map nrb = args(1).map_value();

  if (!error_state) 
    {

      const Cell knots = nrb.contents("knots")(0).cell_value();
      const NDArray coefs = nrb.contents("coefs")(0).array_value();
      octave_idx_type m   = (nrb.contents("number")(0).vector_value())(0) - 1; // m    = size (nrb.coefs, 2) -1;
      octave_idx_type n   = (nrb.contents("number")(0).vector_value())(1) - 1; // n    = size (nrb.coefs, 3) -1;
      octave_idx_type p = (nrb.contents("order")(0).vector_value())(0) - 1;    // p    = nrb.order(1) -1;
      octave_idx_type q = (nrb.contents("order")(0).vector_value())(1) - 1;    // q    = nrb.order(2) -1;

      Array<idx_vector> idx(dim_vector (2, 1), idx_vector(':')); 
      idx(0) = 0;
      const NDArray u(points.index (idx).squeeze ()); // u = points(1,:);

      idx(0) = 1;
      const NDArray v(points.index (idx).squeeze ()); // v = points(2,:);      

      octave_idx_type npt = u.length (); // npt = length(u);
      RowVector M(p+1, 0.0), N (q+1, 0.0);
      Matrix RIkJk(npt, (p+1)*(q+1), 0.0);
      Matrix indIkJk(npt, (p+1)*(q+1), 0.0);
      RowVector denom(npt, 0.0);

      const RowVector U(knots(0).row_vector_value ()); // U = nrb.knots{1};

      const RowVector V(knots(1).row_vector_value ()); // V = nrb.knots{2};
      
      Array<idx_vector> idx2(dim_vector (3, 1), idx_vector(':')); idx2(0) = 3;
      NDArray w (coefs.index (idx2).squeeze ()); // w = squeeze(nrb.coefs(4,:,:));
      
      RowVector spu(u);
      for (octave_idx_type ii(0); ii < npt; ii++)
	{
	  spu(ii) = findspan(m, p, u(ii), U);
	} // spu  =  findspan (m, p, u, U); 

      newargs(3) = U; newargs(2) = p; newargs(1) = u; newargs(0) = spu;
      Matrix Ik = feval (std::string("numbasisfun"), newargs, 1)(0).matrix_value (); // Ik = numbasisfun (spu, u, p, U);

      RowVector spv(v);
      for (octave_idx_type ii(0); ii < v.length(); ii++)
	{
	  spv(ii) = findspan(n, q, v(ii), V);
	} // spv  =  findspan (n, q, v, V);

      newargs(3) = V; newargs(2) = q; newargs(1) = v; newargs(0) = spv;
      Matrix Jk = feval (std::string("numbasisfun"), newargs, 1)(0).matrix_value (); // Jk = numbasisfun (spv, v, q, V);

      Matrix NuIkuk(npt, p+1, 0.0);
      for (octave_idx_type ii(0); ii < npt; ii++)
	{
	  basisfun (int(spu(ii)), u(ii), p, U, M);
	  NuIkuk.insert (M, ii, 0);
	} // NuIkuk = basisfun (spu, u, p, U);

      Matrix NvJkvk(v.length (), q+1, 0.0);
      for (octave_idx_type ii(0); ii < npt; ii++)
	{
	  basisfun(int(spv(ii)), v(ii), q, V, N);
	  NvJkvk.insert (N, ii, 0);
	} // NvJkvk = basisfun (spv, v, q, V);


      for (octave_idx_type k(0); k < npt; k++) 
	for (octave_idx_type ii(0); ii < p+1; ii++) 
	  for (octave_idx_type jj(0); jj < q+1; jj++) 
	    denom(k) += NuIkuk(k, ii) * NvJkvk(k, jj) * w(Ik(k, ii), Jk(k, jj));

      
      for (octave_idx_type k(0); k < npt; k++) 
	for (octave_idx_type ii(0); ii < p+1; ii++) 
	  for (octave_idx_type jj(0); jj < q+1; jj++) 
	    {
	      RIkJk(k, octave_idx_type(ii+(p+1)*jj))  = NuIkuk(k, ii)*NvJkvk(k, jj) * w(Ik(k, ii), Jk(k, jj))/denom(k); 
	      indIkJk(k, octave_idx_type(ii+(p+1)*jj))= Ik(k, ii)+(m+1)*Jk(k, jj)+1;
	    }

      // for k=1:npt
      //       [Jkb, Ika] = meshgrid(Jk(k, :), Ik(k, :)); 
      //       indIkJk(k, :)    = sub2ind([m+1, n+1], Ika(:)+1, Jkb(:)+1);
      //       wIkaJkb(1:p+1, 1:q+1) = reshape (w(indIkJk(k, :)), p+1, q+1); 
      
      //       NuIkukaNvJkvk(1:p+1, 1:q+1) = (NuIkuk(k, :).' * NvJkvk(k, :));
      //       RIkJk(k, :) = (NuIkukaNvJkvk .* wIkaJkb ./ sum(sum(NuIkukaNvJkvk .* wIkaJkb)))(:).';
      //     end
      
      retval(0) = RIkJk; // B = RIkJk;
      retval(1) = indIkJk; // N = indIkJk;

    }
  return retval;
}
