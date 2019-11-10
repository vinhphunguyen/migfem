function ang = vecangle(num,den)

% 
% VECANGLE: An alternative to atan, returning an arctangent in the 
%             range 0 to 2*pi.
% 
% Calling Sequence:
% 
%   ang = vecmag2(num,dum)
% 
% INPUT:
% 
%   num		: Numerator, vector of size (1,nv).
%   dem		: Denominator, vector of size (1,nv).
%
% OUTPUT:
%   ang		: Arctangents, row vector of angles.
% 
% Description:
% 
%   The components of the vector ang are the arctangent of the corresponding
%   enties of num./dem. This function is an alternative for 
%   atan, returning an angle in the range 0 to 2*pi.
% 
% Examples:
% 
%   Find the atan(1.2,2.0) and atan(1.5,3.4) using vecangle
% 
%   ang = vecangle([1.2 1.5], [2.0 3.4]);
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

ang = atan2(num,den);
index = find(ang < 0.0);
ang(index) = 2*pi+ang(index);

end
