function dersv = basisfunder (ii, pl, uu, u_knotl, nders)

% BASISFUNDER:  B-Spline Basis function derivatives.
%
% Calling Sequence:
% 
%   ders = basisfunder (ii, pl, uu, k, nd)
%
%    INPUT:
%   
%      ii  - knot span index (see findspan)
%      pl  - degree of curve
%      uu  - parametric points
%      k   - knot vector
%      nd  - number of derivatives to compute
%
%    OUTPUT:
%   
%      ders - ders(n, i, :) (i-1)-th derivative at n-th point
%   
%    Adapted from Algorithm A2.3 from 'The NURBS BOOK' pg72.
%
%    Copyright (C) 2009 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

  for jj = 1:numel(uu)

    i = ii(jj)+1; %% convert to base-1 numbering of knot spans
    u = uu(jj);

    ders = zeros(nders+1,pl+1);
    ndu = zeros(pl+1,pl+1);
    left = zeros(pl+1);
    right = zeros(pl+1);
    a = zeros(2,pl+1);
    ndu(1,1) = 1;
    for j = 1:pl
      left(j+1) = u - u_knotl(i+1-j);
      right(j+1) = u_knotl(i+j) - u;
      saved = 0;
      for r = 0:j-1
        ndu(j+1,r+1) = right(r+2) + left(j-r+1);
        temp = ndu(r+1,j)/ndu(j+1,r+1);
        ndu(r+1,j+1) = saved + right(r+2)*temp;
        saved = left(j-r+1)*temp;
      end
      ndu(j+1,j+1) = saved;
    end   
    for j = 0:pl
      ders(1,j+1) = ndu(j+1,pl+1);
    end
    for r = 0:pl
      s1 = 0;
      s2 = 1;
      a(1,1) = 1;
      for k = 1:nders %compute kth derivative
        d = 0;
        rk = r-k;
        pk = pl-k;
        if (r >= k)
          a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1);
          d = a(s2+1,1)*ndu(rk+1,pk+1);
        end
        if (rk >= -1)
          j1 = 1;
        else 
          j1 = -rk;
        end
        if ((r-1) <= pk)
          j2 = k-1;
        else 
          j2 = pl-r;
        end
        for j = j1:j2
          a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j))/ndu(pk+2,rk+j+1);
          d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1);
        end
        if (r <= pk)
          a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1);
          d = d + a(s2+1,k+1)*ndu(r+1,pk+1);
        end
        ders(k+1,r+1) = d;
        j = s1;
        s1 = s2;
        s2 = j;
      end
    end
    r = pl;
    for k = 1:nders
      for j = 0:pl
        ders(k+1,j+1) = ders(k+1,j+1)*r;
      end
      r = r*(pl-k);
    end

    dersv(jj, :, :) = ders;
    
  end

end

%!test
%! k    = [0 0 0 0 1 1 1 1];
%! p    = 3;
%! u    = rand (1);
%! i    = findspan (numel(k)-p-2, p, u, k);
%! ders = basisfunder (i, p, u, k, 1);
%! sumders = sum (squeeze(ders), 2);
%! assert (sumders(1), 1, 1e-15);
%! assert (sumders(2:end), 0, 1e-15);

%!test
%! k    = [0 0 0 0 1/3 2/3 1 1 1 1];
%! p    = 3;
%! u    = rand (1);
%! i    = findspan (numel(k)-p-2, p, u, k);
%! ders = basisfunder (i, p, u, k, 7); 
%! sumders = sum (squeeze(ders), 2);
%! assert (sumders(1), 1, 1e-15);
%! assert (sumders(2:end), zeros(rows(squeeze(ders))-1, 1), 1e-13);

%!test
%! k    = [0 0 0 0 1/3 2/3 1 1 1 1];
%! p    = 3;
%! u    = rand (100, 1);
%! i    = findspan (numel(k)-p-2, p, u, k);
%! ders = basisfunder (i, p, u, k, 7);
%! for ii=1:10
%!   sumders = sum (squeeze(ders(ii,:,:)), 2);
%!   assert (sumders(1), 1, 1e-15);
%!   assert (sumders(2:end), zeros(rows(squeeze(ders(ii,:,:)))-1, 1), 1e-13);
%! end
%! assert (ders(:, (p+2):end, :), zeros(numel(u), 8-p-1, p+1), 1e-13)
%! assert (all(all(ders(:, 1, :) <= 1)), true)



