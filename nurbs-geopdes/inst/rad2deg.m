function deg = rad2deg(rad)
% 
% RAD2DEG: Convert radians to degrees.
% 
% Calling Sequence:
% 
%   rad = rad2deg(deg);
% 
% INPUT:
% 
%   rad		: Angle in radians.
%
% OUTPUT:
%
%   deg		: Angle in degrees.
% 
% Description:
% 
%   Convenient utility function for converting radians to degrees, which are
%   often the required angular units for functions in the NURBS toolbox.
% 
% Examples:
% 
%   Convert 0.3 radians to degrees
% 
%   rad = deg2rad(0.3);
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

deg = 180.0*rad/pi;

end
