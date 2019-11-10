function ck = curvederiveval (n, p, U, P, u, d) 
%
% CURVEDERIVEVAL: Compute the derivatives of a B-spline curve.
% 
% usage: ck = curvederiveval (n, p, U, P, u, d) 
%
%  INPUT: 
%
%        n+1 = number of control points
%        p   = spline order
%        U   = knots
%        P   = control points
%        u   = evaluation point
%        d   = derivative order
%
%  OUTPUT:
%
%        ck (k+1) =  curve differentiated k times
%
% Adaptation of algorithm A3.4 from the NURBS book, pg99
%
%    Copyright (C) 2009 Carlo de Falco
%    Copyright (C) 2010 Rafael Vazquez
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

  ck = zeros (d+1, 1);
  du = min (d, p);   

  span = findspan (n, p, u, U);
  for ip=0:p
      N(1:ip+1,ip+1) = basisfun (span, u, ip, U)';
  end
    
  pk = curvederivcpts (n, p, U, P, du, span-p, span);

  for k = 0:du
    for j = 0:p-k
      ck(k+1) = ck(k+1) + N(j+1,p-k+1)*pk(k+1,j+1);
    end
  end
  
end

%!test
%! k = [0 0 0 1 1 1];
%! coefs(:,1) = [0;0;0;1];
%! coefs(:,2) = [1;0;1;1];
%! coefs(:,3) = [1;1;1;1];
%! crv = nrbmak (coefs, k);
%! ck = curvederiveval (crv.number-1, crv.order-1, crv.knots, squeeze (crv.coefs(1,:,:)), 0.5, 2);
%! assert(ck, [0.75; 1; -2]);
%! ck = curvederiveval (crv.number-1, crv.order-1, crv.knots, squeeze (crv.coefs(2,:,:)), 0.5, 2);
%! assert(ck, [0.25; 1; 2]);
%! ck = curvederiveval (crv.number-1, crv.order-1, crv.knots, squeeze (crv.coefs(3,:,:)), 0.5, 2);
%! assert(ck, [0.75; 1; -2]);
