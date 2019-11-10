function ck = nrbcrvderiveval (crv, u, d) 

%
% NRBCRVDERIVEVAL: Evaluate n-th order derivatives of a NURBS curve.
%
% usage: skl = nrbcrvderiveval (crv, u, d) 
%
%   INPUT:
%
%   crv : NURBS curve structure, see nrbmak
%
%   u   : parametric coordinate of the points where we compute the derivatives
%
%   d   : number of partial derivatives to compute
%
%
%   OUTPUT: 
%
%   ck (i, j, l) = i-th component derived j-1 times at the l-th point.
%
% Adaptation of algorithm A4.2 from the NURBS book, pg127
%
%    Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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
  
 ck = zeros (3, d+1, numel(u));
 
 for iu = 1:numel(u);
   wders = squeeze (curvederiveval (crv.number-1, crv.order-1,  ...
				     crv.knots, squeeze (crv.coefs(4, :)), u(iu), d));

   for idim = 1:3
     Aders = squeeze (curvederiveval (crv.number-1, crv.order-1,  ...
			crv.knots, squeeze (crv.coefs(idim, :)), u(iu), d));
     for k=0:d
       v = Aders(k+1);
       for i=1:k
	     v = v - nchoosek(k,i)*wders(i+1)*ck(idim, k-i+1, iu);
       end
	   ck(idim, k+1, iu) = v/wders(1);
     end
   end
 end
end

%!test
%! knots = [0 0 0 1 1 1];
%! coefs(:,1) = [0; 0; 0; 1];
%! coefs(:,2) = [1; 0; 1; 1];
%! coefs(:,3) = [1; 1; 1; 2];
%! crv = nrbmak (coefs, knots);
%! u = linspace (0, 1, 10);
%! ck = nrbcrvderiveval (crv, u, 2);
%! w  = @(x) 1 + x.^2;
%! dw = @(x) 2*x;
%! F1 = @(x) (2*x - x.^2)./w(x);
%! F2 = @(x) x.^2./w(x);
%! F3 = @(x) (2*x - x.^2)./w(x);
%! dF1 = @(x) (2 - 2*x)./w(x) - 2*(2*x - x.^2).*x./w(x).^2;
%! dF2 = @(x) 2*x./w(x) - 2*x.^3./w(x).^2;
%! dF3 = @(x) (2 - 2*x)./w(x) - 2*(2*x - x.^2).*x./w(x).^2;
%! d2F1 = @(x) -2./w(x) - 2*x.*(2-2*x)./w(x).^2 - (8*x-6*x.^2)./w(x).^2 + 8*x.^2.*(2*x-x.^2)./w(x).^3;
%! d2F2 = @(x) 2./w(x) - 4*x.^2./w(x).^2 - 6*x.^2./w(x).^2 + 8*x.^4./w(x).^3;
%! d2F3 = @(x) -2./w(x) - 2*x.*(2-2*x)./w(x).^2 - (8*x-6*x.^2)./w(x).^2 + 8*x.^2.*(2*x-x.^2)./w(x).^3;
%! assert ([F1(u); F2(u); F3(u)], squeeze(ck(:, 1, :)), 1e2*eps);
%! assert ([dF1(u); dF2(u); dF3(u)], squeeze(ck(:, 2, :)), 1e2*eps);
%! assert ([d2F1(u); d2F2(u); d2F3(u)], squeeze(ck(:, 3, :)), 1e2*eps);
