function s = findspan(n,p,u,U)                
% FINDSPAN  Find the span of a B-Spline knot vector at a parametric point
%
% Calling Sequence:
% 
%   s = findspan(n,p,u,U)
% 
%  INPUT:
% 
%    n - number of control points - 1
%    p - spline degree
%    u - parametric point
%    U - knot sequence
% 
%  OUTPUT:
% 
%    s - knot span index
%
%  Modification of Algorithm A2.1 from 'The NURBS BOOK' pg68
%
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

if (max(u(:))>U(end) || min(u(:))<U(1))
  error('Some value is outside the knot span')
end

s = zeros(size(u));
for j = 1:numel(u)
  if (u(j)==U(n+2)), s(j)=n; continue, end
  s(j) = find(u(j) >= U,1,'last')-1;
end

end

%!test
%!  n = 3; 
%!  U = [0 0 0 1/2 1 1 1]; 
%!  p = 2; 
%!  u = linspace(0, 1, 10);  
%!  s = findspan (n, p, u, U);
%!  assert (s, [2*ones(1, 5) 3*ones(1, 5)]);

%!test
%! p = 2; m = 7; n = m - p - 1;
%! U = [zeros(1,p)  linspace(0,1,m+1-2*p) ones(1,p)];
%! u = [ 0   0.11880   0.55118   0.93141   0.40068   0.35492 0.44392   0.88360   0.35414   0.92186   0.83085   1];
%! s = [2   2   3   4   3   3   3   4   3   4   4   4];
%! assert (findspan (n, p, u, U), s, 1e-10);
