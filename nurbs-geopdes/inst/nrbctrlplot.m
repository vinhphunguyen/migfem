function nrbctrlplot (nurbs)

% NRBCTRLPLOT: Plot a NURBS entity along with its control points.
% 
% Calling Sequence:
% 
%   nrbctrlplot (nurbs)
% 
% INPUT:
% 
%   nurbs: NURBS curve or surface, see nrbmak.
% 
% Example:
%
%   Plot the test curve and test surface with their control polygon and
%    control net, respectively
%
%   nrbctrlplot(nrbtestcrv)
%   nrbctrlplot(nrbtestsrf)
%
% See also:
% 
%   nrbkntplot
%
%    Copyright (C) 2011 Rafael Vazquez
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

if (nargin < 1)
  error ('nrbctrlplot: Need a NURBS to plot!');
end

% Default values
light='on';
cmap='summer';

colormap (cmap);

hold_flag = ishold;

if (iscell (nurbs.knots))
  if (size (nurbs.knots,2) == 3)
    error ('nrbctrlplot: not implemented for NURBS volumes');
  elseif (size (nurbs.knots,2) == 2) % plot a NURBS surface

    nsub = 100;
    nrbplot (nurbs, [nsub nsub], 'light', light, 'colormap', cmap);
    hold on

% And plot the control net
    for ii = 1:size (nurbs.coefs, 2)
      coefs = reshape (nurbs.coefs(1:3,ii,:), 3, []);
      weights = reshape (nurbs.coefs(4,ii,:), 1, []);
      plot3 (coefs(1,:)./weights, coefs(2,:)./weights, coefs(3,:)./weights,'k--')
      plot3 (coefs(1,:)./weights, coefs(2,:)./weights, coefs(3,:)./weights,'r.','MarkerSize',20)
    end
    for jj = 1:size (nurbs.coefs, 3)
      coefs = reshape (nurbs.coefs(1:3,:,jj), 3, []);
      weights = reshape (nurbs.coefs(4,:,jj), 1, []);
      plot3 (coefs(1,:)./weights, coefs(2,:)./weights, coefs(3,:)./weights,'k--')
      plot3 (coefs(1,:)./weights, coefs(2,:)./weights, coefs(3,:)./weights,'r.','MarkerSize',20)
    end
  end
else % plot a NURBS curve
  nsub = 1000;
  nrbplot (nurbs, nsub);
  hold on

% And plot the control polygon
  coefs = nurbs.coefs(1:3,:);
  weights = nurbs.coefs(4,:);
  plot3 (coefs(1,:)./weights, coefs(2,:)./weights, coefs(3,:)./weights,'k--')
  plot3 (coefs(1,:)./weights, coefs(2,:)./weights, coefs(3,:)./weights,'r.','MarkerSize',20)
end

if (~hold_flag)
  hold off
end

end

%!demo
%! crv = nrbtestcrv;
%! nrbctrlplot(crv)
%! title('Test curve')
%! hold off

%!demo
%! srf = nrbtestsrf;
%! nrbctrlplot(srf)
%! title('Test surface')
%! hold off

