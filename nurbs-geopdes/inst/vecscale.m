function ss = vecscale(vector)

%
% VECSCALE: Transformation matrix for a scaling.
% 
% Calling Sequence:
% 
%   ss = vecscale(svec)
% 
% INPUT:
% 
%   svec    : A vectors defining the scaling along the x,y and z axes.
%             i.e. [sx, sy, sy]
%
% OUTPUT:
% 
%   ss	    : Scaling Transformation Matrix
% 
% Description:
% 
%   Returns a (4x4) Transformation matrix for scaling.
% 
%   The matrix is:
% 
%         [ sx  0   0   0]
%         [ 0   sy  0   0]
%         [ 0   0   sz  0]
%         [ 0   0   0   1]
% 
% Example:
% 
%   Scale up the NURBS line (0.0,0.0,0.0) - (1.0,1.0,1.0) by 3 along
%   the x-axis, 2 along the y-axis and 4 along the z-axis.
% 
%   line = nrbline([0.0 0.0 0.0],[1.0 1.0 1.0]);
%   trans = vecscale([3.0 2.0 4.0]);
%   sline = nrbtform(line, trans);
% 
% See also:
% 
%    nrbtform
%
%    Copyright (C) 2000 Mark Spink, 2007 Daniel Claxton
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
  error('Scaling vector not specified');
end   

s = [vector(:);0;0];
ss = [s(1) 0 0 0; 0 s(2) 0 0; 0 0 s(3) 0; 0 0 0 1];

end
