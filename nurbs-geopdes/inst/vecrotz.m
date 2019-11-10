function rz = vecrotz(angle)
% 
% VECROTZ: Transformation matrix for a rotation around the z axis. 
% 
% Calling Sequence:
% 
%   rz = vecrotz(angle);
% 
% INPUT:
% 
%   angle	: rotation angle defined in radians
%
% OUTPUT:
% 
%   rz		: (4x4) Transformation matrix.
% 
% 
% Description:
% 
%   Return the (4x4) Transformation matrix for a rotation about the z axis
%   by the defined angle.
% 
%   The matrix is:
% 
%         [  cos(angle)  -sin(angle)       0          0]
%         [ -sin(angle)   cos(angle)       0          0]
%         [      0            0            1          0]
%         [      0            0            0          1]
% 
% Examples:
% 
%  Rotate the NURBS line (0.0 0.0 0.0) - (3.0 3.0 3.0) by 45 degrees
%  around the z-axis
% 
%    line = nrbline([0.0 0.0 0.0],[3.0 3.0 3.0]);
%    trans = vecrotz(%pi/4);
%    rline = nrbtform(line, trans);
% 
% See also:
% 
%    nrbtform
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

sn = sin(angle);
cn = cos(angle);
rz = [cn -sn 0 0; sn cn 0 0; 0 0 1 0; 0 0 0 1];

end
