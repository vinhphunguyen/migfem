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

octave_idx_type findspan(int n, int p, double u, const RowVector& U);

void basisfun(int i, double u, int p, const RowVector& U, RowVector& N);

void basisfunder (int i, int pl, double uu, const RowVector& u_knotl, 
		  int nders, NDArray& dersv);

int curvederivcpts (octave_idx_type n, octave_idx_type p, 
		    const RowVector &U, const NDArray &P, 
		    octave_idx_type  d, octave_idx_type r1, 
		    octave_idx_type r2, 
		    Matrix &pk);

int surfderivcpts (octave_idx_type n, octave_idx_type  p, const RowVector& U, 
		   octave_idx_type m, octave_idx_type q, const RowVector& V, 
		   const Matrix& P, octave_idx_type d, octave_idx_type r1, 
		   octave_idx_type r2, octave_idx_type s1, 
		   octave_idx_type s2, NDArray &pkl);

int surfderiveval (octave_idx_type n, octave_idx_type p, const RowVector &U, 
		   octave_idx_type m, octave_idx_type q, const RowVector &V, 
		   const Matrix &P, double u, double v, octave_idx_type d, 
		   Matrix &skl);
