function N = tbasisfun (u, p, U)
%
% TBASISFUN: Compute a B- or T-Spline basis function from its local knot vector.
%
% usage:
%
% N = tbasisfun (u, p, U)
% N = tbasisfun ([u; v], [p q], {U, V})
% N = tbasisfun ([u; v; w], [p q r], {U, V, W})
% 
% INPUT:
%
%  u or [u; v] : points in parameter space where the basis function is to be
%  evaluated 
%  
%  U or {U, V} : local knot vector
%
% p or [p q] : polynomial order of the basis function
%
% OUTPUT:
%
%  N : basis function evaluated at the given parametric points 
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
  
  if (~ iscell (U))
    U = sort (U);
    assert (numel (U) == p+2)
    
    N = zeros(1,numel(u));
    for ii=1:numel(u)
      N(ii) = onebasisfun__ (u(ii), p, U);
    end

  elseif size(U,2) == 2
    U{1} = sort(U{1}); U{2} = sort(U{2});
    assert (numel(U{1}) == p(1)+2 && numel(U{2}) == p(2)+2)
    
    Nu = zeros(1,numel(u(1,:))); Nv = zeros(1,numel(u(1,:)));
    for ii=1:numel(u(1,:))
      Nu(ii) = onebasisfun__ (u(1,ii), p(1), U{1});
    end

    for ii=1:numel(u(1,:))
      Nv(ii) = onebasisfun__ (u(2,ii), p(2), U{2});
    end

    N = Nu.*Nv;

  elseif size(U,2) == 3
    U{1} = sort(U{1}); U{2} = sort(U{2}); U{3} = sort(U{3});
    assert (numel(U{1}) == p(1)+2 && numel(U{2}) == p(2)+2 && numel(U{3}) == p(3)+2)
    
    Nu = zeros(1,numel(u(1,:))); Nv = zeros(1,numel(u(1,:))); Nw = zeros(1,numel(u(1,:)));
    for ii=1:numel(u(1,:))
      Nu(ii) = onebasisfun__ (u(1,ii), p(1), U{1});
    end

    for ii=1:numel(u(1,:))
      Nv(ii) = onebasisfun__ (u(2,ii), p(2), U{2});
    end

    for ii=1:numel(u(1,:))
      Nw(ii) = onebasisfun__ (u(3,ii), p(3), U{3});
    end

    N = Nu.*Nv.*Nw;
  end

end

%!demo
%! U = {[0 0 1/2 1 1], [0 0 0 1 1]};
%! p = [3, 3];
%! [X, Y] = meshgrid (linspace(0, 1, 30));
%! u = [X(:), Y(:)]';
%! N = tbasisfun (u, p, U);
%! surf (X, Y, reshape (N, size(X)))
%! title('Basis function associated to a local knot vector')
%! hold off
