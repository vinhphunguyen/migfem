function srf = nrbtestsrf
% NRBTESTSRF: Constructs a simple test surface.
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

% allocate multi-dimensional array of control points
pnts = zeros(3,5,5);

% define a grid of control points
% in this case a regular grid of u,v points
% pnts(3,u,v)
%

pnts(:,:,1) = [ 0.0  3.0  5.0  8.0 10.0;     % w*x
                0.0  0.0  0.0  0.0  0.0;     % w*y
                2.0  2.0  7.0  7.0  8.0];    % w*z

pnts(:,:,2) = [ 0.0  3.0  5.0  8.0 10.0;
                3.0  3.0  3.0  3.0  3.0;
                0.0  0.0  5.0  5.0  7.0];

pnts(:,:,3) = [ 0.0  3.0  5.0  8.0 10.0;
                5.0  5.0  5.0  5.0  5.0;
                0.0  0.0  5.0  5.0  7.0];

pnts(:,:,4) = [ 0.0  3.0  5.0  8.0 10.0;
                8.0  8.0  8.0  8.0  8.0;
                5.0  5.0  8.0  8.0 10.0];

pnts(:,:,5) = [ 0.0  3.0  5.0  8.0 10.0;
               10.0 10.0 10.0 10.0 10.0;
                5.0  5.0  8.0  8.0 10.0];

% knots
knots{1} = [0 0 0 1/3 2/3 1 1 1]; % knots along u
knots{2} = [0 0 0 1/3 2/3 1 1 1]; % knots along v

% make and draw nurbs surface
srf = nrbmak(pnts,knots);

end

%!demo
%! srf = nrbtestsrf;
%! nrbplot(srf,[20 30])
%! title('Test surface')
%! hold off

