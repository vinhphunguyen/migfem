function curve = nrbrect(w,h)
% 
% NRBRECT: Construct NURBS representation of a rectangular curve.
% 
% Calling Sequence:
% 
%   crv = nrbrect()
%   crv = nrbrect(size)
%   crv = nrbrect(width, height)
% 
% INPUT:
% 
%   size	: Size of the square (width = height).
% 
%   width	: Width of the rectangle (along x-axis).
% 
%   height	: Height of the rectangle (along y-axis).
%
% OUTPUT:
%
%   crv		: NURBS curve, see nrbmak. 
%  
% 
% Description:
% 
%   Construct a rectangle or square in the x-y plane with the bottom
%   lhs corner at (0,0,0). If no rhs arguments provided the function
%   constructs a unit square.
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

if nargin < 1
   w = 1;
   h = 1;
end

if nargin < 2
   h = w;
end

coefs  = [0 w w w w 0 0 0;
          0 0 0 h h h h 0;
          0 0 0 0 0 0 0 0;
          1 1 1 1 1 1 1 1];

knots  = [0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1];

curve = nrbmak(coefs, knots);

end

%!demo
%! crv = nrbtform(nrbrect(2,1), vecrotz(deg2rad(35)));
%! nrbplot(crv,4);
%! axis equal
%! title('Construction and rotation of a rectangular curve.');
%! hold off