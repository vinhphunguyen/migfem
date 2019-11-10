function dd = vectrans(vector)
% 
% VECTRANS: Transformation matrix for a translation.
% 
% Calling Sequence:
% 
%   st = vectrans(tvec)
% 
% INPUT:
% 
%   tvec	: A vectors defining the translation along the x,y and
%                   z axes. i.e. [tx, ty, ty]
%
% OUTPUT:
% 
%   st		: Translation Transformation Matrix
% 
% Description:
% 
%   Returns a (4x4) Transformation matrix for translation.
% 
%   The matrix is:
% 
%         [ 1   0   0   tx ]
%         [ 0   0   0   ty ]
%         [ 0   0   0   tz ]
%         [ 0   0   0   1  ]
% 
% Examples:
% 
%   Translate the NURBS line (0.0,0.0,0.0) - (1.0,1.0,1.0) by 3 along
%   the x-axis, 2 along the y-axis and 4 along the z-axis.
%
%   line = nrbline([0.0 0.0 0.0],[1.0 1.0 1.0]);
%   trans = vectrans([3.0 2.0 4.0]);
%   tline = nrbtform(line, trans);
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


if nargin < 1
  error('Translation vector required');
end   

v = [vector(:);0;0];
dd = [1 0 0 v(1); 0 1 0 v(2); 0 0 1 v(3); 0 0 0 1];

end
