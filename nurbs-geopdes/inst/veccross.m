function cross = veccross(vec1,vec2)
% 
% VECCROSS: The cross product of two vectors.
% 
% Calling Sequence:
% 
%   cross = veccross(vec1,vec2);
% 
% INPUT:
% 
%   vec1	: An array of column vectors represented by a matrix of
%   vec2	size (dim,nv), where is the dimension of the vector and
% 		nv the number of vectors.
%
% OUTPUT:
% 
%   cross	: Array of column vectors, each element is corresponding
% 		to the cross product of the respective components in vec1
% 		and vec2.
% 
% Description:
% 
%   Cross product of two vectors.
% 
% Examples:
% 
%   Determine the cross products of:
%   (2.3,3.4,5.6) and (1.2,4.5,1.2)
%   (5.1,0.0,2.3) and (2.5,3.2,4.0)
% 
%   cross = veccross([2.3 5.1; 3.4 0.0; 5.6 2.3],[1.2 2.5; 4.5 3.2; 1.2 4.0]);
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

if size(vec1,1) == 2
  % 2D vector
  cross = zeros(size(vec1));
  cross(3,:) = vec1(1,:).*vec2(2,:)-vec1(2,:).*vec2(1,:);
else
  % 3D vector
  cross = [vec1(2,:).*vec2(3,:)-vec1(3,:).*vec2(2,:);
           vec1(3,:).*vec2(1,:)-vec1(1,:).*vec2(3,:);
           vec1(1,:).*vec2(2,:)-vec1(2,:).*vec2(1,:)];
end

end
