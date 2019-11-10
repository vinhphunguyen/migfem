function curve = nrbline(p1,p2)
% 
% NRBLINE: Construct a straight line.
% 
% Calling Sequence:
% 
%   crv = nrbline()
%   crv = nrbline(p1,p2)
% 
% INPUT:
% 
% p1		: 2D or 3D cartesian coordinate of the start point.
% 
% p2            : 2D or 3D cartesian coordinate of the end point.
%
% OUTPUT:
% 
% crv		: NURBS curve for a straight line.
% 
% Description:
% 
%   Constructs NURBS data structure for a straight line. If no rhs 
%   coordinates are included the function returns a unit straight
%   line along the x-axis.
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

coefs = [zeros(3,2); ones(1,2)];

if nargin < 2
  coefs(1,2) = 1.0;  
else
  coefs(1:length(p1),1) = p1(:);    
  coefs(1:length(p2),2) = p2(:);
end

curve = nrbmak(coefs, [0 0 1 1]);

end

%!demo
%! crv = nrbline([0.0 0.0 0.0]',[5.0 4.0 2.0]');
%! nrbplot(crv,1);
%! grid on;
%! title('3D straight line.');
%! hold off
