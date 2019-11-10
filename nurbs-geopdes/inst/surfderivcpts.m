function pkl = surfderivcpts (n, p, U, m, q, V, P, d, r1, r2, s1, s2) 
%
% SURFDERIVCPTS: Compute control points of n-th derivatives of a NURBS surface.
% 
% usage: pkl = surfderivcpts (n, p, U, m, q, V, P, d) 
%
%  INPUT: 
%
%        n+1, m+1 = number of control points
%        p, q     = spline order
%        U, V     = knots
%        P        = control points
%        d        = derivative order
%
%  OUTPUT:
%
%        pkl (k+1, l+1, i+1, j+1) = i,jth control point
%                                   of the surface differentiated k
%                                   times in the u direction and l
%                                   times in the v direction
%
% Adaptation of algorithm A3.7 from the NURBS book, pg114
%
%    Copyright (C) 2009 Carlo de Falco
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
  
  if (nargin <= 8)
    r1 = 0; r2 = n;
    s1 = 0; s2 = m;
  end
  r = r2-r1;
  s = s2-s1;

  du = min (d, p);   dv = min (d, q); 

  for j=s1:s2
    temp = curvederivcpts (n, p, U, P(:,j+1:end), du, r1, r2);
    for k=0:du
      for i=0:r-k
       pkl (k+1, 1, i+1, j-s1+1) = temp (k+1, i+1);
      end
    end
  end
  
  for k=0:du
    for i=0:r-k
      dd = min (d-k, dv);
      temp = curvederivcpts (m, q, V(s1+1:end), pkl(k+1, 1, i+1, :),  ...
			     dd, 0, s);
      for l=1:dd
       for j=0:s-l
        pkl (k+1, l+1, i+1, j+1) = temp (l+1, j+1);
       end
      end
    end
  end

end

%!test
%! coefs = cat(3,[0 0; 0 1],[1 1; 0 1]);
%! knots = {[0 0 1 1]  [0 0 1 1]};
%! plane = nrbmak(coefs,knots);
%! pkl = surfderivcpts (plane.number(1)-1, plane.order(1)-1,...
%!                       plane.knots{1}, plane.number(2)-1,...
%!                       plane.order(2)-1, plane.knots{2}, ...
%!                       squeeze (plane.coefs(1,:,:)), 1);
