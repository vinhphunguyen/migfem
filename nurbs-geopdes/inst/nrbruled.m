function srf = nrbruled (crv1, crv2)

% NRBRULED: Construct a ruled surface between two NURBS curves.
% 
% Calling Sequence:
% 
%   srf = nrbruled(crv1, crv2)
% 
% INPUT:
% 
%   crv1	: First NURBS curve, see nrbmak.
% 
%   crv2	: Second NURBS curve, see nrbmak.
%
% OUTPUT:
% 
%   srf		: Ruled NURBS surface.
% 
% Description:
% 
%   Constructs a ruled surface between two NURBS curves. The ruled surface is
%   ruled along the V direction.
% 
% Examples:
% 
%   Construct a ruled surface between a semicircle and a straight line.
% 
%   cir = nrbcirc(1,[0 0 0],0,pi);
%   line = nrbline([-1 0.5 1],[1 0.5 1]);
%   srf = nrbruled(cir,line);
%   nrbplot(srf,[20 20]);
%
%    Copyright (C) 2000 Mark Spink
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

if (iscell(crv1.knots) || iscell(crv2.knots))
  error ('Both NURBS must be curves');
end

% ensure both curves have a common degree
d = max ([crv1.order, crv2.order]);
crv1 = nrbdegelev (crv1, d - crv1.order);
crv2 = nrbdegelev (crv2, d - crv2.order);

% merge the knot vectors, to obtain a common knot vector
k1 = crv1.knots;
k2 = crv2.knots;
ku = unique ([k1 k2]);
n = length (ku);
ka = [];
kb = [];
for i = 1:n
  i1 = length (find (k1 == ku(i)));
  i2 = length (find (k2 == ku(i)));
  m = max (i1, i2);
  ka = [ka ku(i)*ones(1,m-i1)];
  kb = [kb ku(i)*ones(1,m-i2)];
end
crv1 = nrbkntins (crv1, ka);
crv2 = nrbkntins (crv2, kb);

coefs(:,:,1) = crv1.coefs;
coefs(:,:,2) = crv2.coefs;
srf = nrbmak (coefs, {crv1.knots, [0 0 1 1]});

end

%!demo
%! pnts = [0.5 1.5 4.5 3.0 7.5 6.0 8.5;
%!         3.0 5.5 5.5 1.5 1.5 4.0 4.5;
%!         0.0 0.0 0.0 0.0 0.0 0.0 0.0];
%! crv1 = nrbmak (pnts,[0 0 0 1/4 1/2 3/4 3/4 1 1 1]);
%! crv2 = nrbtform (nrbcirc (4,[4.5;0],pi,0.0),vectrans([0.0 4.0 -4.0]));
%! srf = nrbruled (crv1,crv2);
%! nrbplot (srf,[40 20]);
%! title ('Ruled surface construction from two NURBS curves.');
%! hold off