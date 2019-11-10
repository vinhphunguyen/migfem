function B = numbasisfun (iv, uv, p, U)

% NUMBASISFUN:  List non-zero Basis functions for B-Spline in a given knot-span
%
% Calling Sequence:
% 
%   N = numbasisfun(i,u,p,U)
%   
%    INPUT:
%   
%      i - knot span  ( from FindSpan() )
%      u - parametric point
%      p - spline degree
%      U - knot sequence
%   
%    OUTPUT:
%   
%      N - Basis functions (numel(u)x(p+1))
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

B = bsxfun (@(a, b) a+b,iv-p, (0:p).').';

end

%!test
%!  n = 3; 
%!  U = [0 0 0 1/2 1 1 1]; 
%!  p = 2; 
%!  u = linspace (0, 1, 10);  
%!  s = findspan (n, p, u, U); 
%!  Bref = [0   0   0   0   0   1   1   1   1   1; ...
%!          1   1   1   1   1   2   2   2   2   2; ...
%!          2   2   2   2   2   3   3   3   3   3].';
%!  B = numbasisfun (s, u, p, U);
%!  assert (B, Bref)