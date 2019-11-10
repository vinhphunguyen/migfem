function surf = nrbcylind(height,radius,center,sang,eang)
% 
% NRBCYLIND: Construct a cylinder or cylindrical patch.
% 
% Calling Sequence:
% 
%   srf = nrbcylind()
%   srf = nrbcylind(height)
%   srf = nrbcylind(height,radius)
%   srf = nrbcylind(height,radius,center)
%   srf = nrbcylind(height,radius,center,sang,eang)
% 
% INPUT:
% 
%   height	: Height of the cylinder along the axis, default 1.0
% 
%   radius	: Radius of the cylinder, default 1.0
% 
%   center	: Center of the cylinder, default (0,0,0)
% 
%   sang	: Start angle relative to the origin, default 0.
% 
%   eang	: End angle relative to the origin, default 2*pi.
%
% OUTPUT: 
%
%   srf     : cylindrical surface patch 
% 
% Description:
% 
%   Construct a cylinder or cylindrical patch by extruding a circular arc.
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
  height = 1;
end

if nargin < 2
  radius = 1;
end

if nargin < 3
  center = [];
end
   
if nargin < 5
  sang = 0;
  eang = 2*pi;
end

surf = nrbextrude(nrbcirc(radius,center,sang,eang),[0.0 0.0 height]);

end

%!demo
%! srf = nrbcylind(3,1,[],deg2rad(270),deg2rad(180));
%! nrbplot(srf,[20,20]);
%! axis equal;
%! title('Cylinderical section by extrusion of a circular arc.');
%! hold off