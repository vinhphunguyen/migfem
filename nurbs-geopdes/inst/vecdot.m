function dot = vecdot(vec1,vec2)
% 
% VECDOT: The dot product of two vectors.
% 
% Calling Sequence:
% 
%   dot = vecdot(vec1,vec2);
% 
% INPUT:
% 
%   vec1	: An array of column vectors represented by a matrix of
%   vec2	size (dim,nv), where is the dimension of the vector and
% 		nv the number of vectors.
%
% OUTPUT:
%
%   dot		: Row vector of scalars, each element corresponding to
% 		the dot product of the respective components in vec1 and
% 		vec2.
% 
% Description:
% 
%   Scalar dot product of two vectors.
% 
% Examples:
% 
%   Determine the dot product of
%   (2.3,3.4,5.6) and (1.2,4.5,1.2)
%   (5.1,0.0,2.3) and (2.5,3.2,4.0)
%
%   dot = vecdot([2.3 5.1; 3.4 0.0; 5.6 2.3],[1.2 2.5; 4.5 3.2; 1.2 4.0]);
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

dot = sum(vec1.*vec2);

end
