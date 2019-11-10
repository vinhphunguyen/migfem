function nrbkntplot (nurbs)

% NRBKNTPLOT: Plot a NURBS entity with the knots subdivision.
% 
% Calling Sequence:
% 
%   nrbkntplot(nurbs)
% 
% INPUT:
% 
%   nurbs: NURBS curve, surface or volume, see nrbmak.
% 
% Example:
%
%   Plot the test surface with its knot vector
%
%   nrbkntplot(nrbtestsrf)
%
% See also:
% 
%   nrbctrlplot
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
  error ('nrbkntplot: Need a NURBS to plot!');
end

% Default values
light='on';
cmap='summer';

colormap (cmap);

hold_flag = ishold;

if (iscell (nurbs.knots))
 if (size (nurbs.knots,2) == 2) % plot a NURBS surface
   nsub = 100;
   nrbplot (nurbs, [nsub nsub], 'light', light, 'colormap', cmap);
   hold on

   % And plot the knots
   knt1 = unique (nurbs.knots{1});
   knt2 = unique (nurbs.knots{2});
   p1 = nrbeval (nurbs, {knt1, linspace(0.0,1.0,nsub)});
   p2 = nrbeval (nurbs, {linspace(0.0,1.0,nsub), knt2});

  if (any (nurbs.coefs(3,:)))
    % surface in a 3D space
    for ii = 1:numel(knt1)
      plot3 (squeeze(p1(1,ii,:)), squeeze(p1(2,ii,:)), squeeze(p1(3,ii,:)));
    end
    for ii = 1:numel(knt2)
      plot3 (squeeze(p2(1,:,ii)), squeeze(p2(2,:,ii)), squeeze(p2(3,:,ii))); 
    end
  else
    % plain surface
    for ii = 1:numel(knt1)
      plot (squeeze(p1(1,ii,:)), squeeze (p1(2,ii,:))); 
    end
    for ii = 1:numel(knt2)
      plot (p2(1,:,ii),p2(2,:,ii));
    end
  end



 elseif (size (nurbs.knots,2) == 3) % plot a NURBS volume
   nsub = 30;
   nrbplot (nurbs, [nsub nsub nsub], 'light', light, 'colormap', cmap);
   hold on

   % And plot the knots
   knt1 = unique (nurbs.knots{1});
   knt2 = unique (nurbs.knots{2});
   knt3 = unique (nurbs.knots{3});
   kv_face1 = nrbeval (nurbs, {0, knt2, linspace(0.0,1.0,nsub)});
   kw_face1 = nrbeval (nurbs, {0, linspace(0.0,1.0,nsub), knt3});
   kv_face2 = nrbeval (nurbs, {1, knt2, linspace(0.0,1.0,nsub)});
   kw_face2 = nrbeval (nurbs, {1, linspace(0.0,1.0,nsub), knt3});
   ku_face3 = nrbeval (nurbs, {knt1, 0, linspace(0.0,1.0,nsub)});
   kw_face3 = nrbeval (nurbs, {linspace(0.0,1.0,nsub), 0, knt3});
   ku_face4 = nrbeval (nurbs, {knt1, 1, linspace(0.0,1.0,nsub)});
   kw_face4 = nrbeval (nurbs, {linspace(0.0,1.0,nsub), 1, knt3});
   ku_face5 = nrbeval (nurbs, {knt1, linspace(0.0,1.0,nsub), 0});
   kv_face5 = nrbeval (nurbs, {linspace(0.0,1.0,nsub), knt2, 0});
   ku_face6 = nrbeval (nurbs, {knt1, linspace(0.0,1.0,nsub), 1});
   kv_face6 = nrbeval (nurbs, {linspace(0.0,1.0,nsub), knt2, 1});

   for ii = 1:numel(knt1)
     plot3 (squeeze (ku_face3(1,ii,:,:)), squeeze (ku_face3(2,ii,:,:)), squeeze (ku_face3(3,ii,:,:))); 
     plot3 (squeeze (ku_face4(1,ii,:,:)), squeeze (ku_face4(2,ii,:,:)), squeeze (ku_face4(3,ii,:,:))); 
     plot3 (squeeze (ku_face5(1,ii,:,:)), squeeze (ku_face5(2,ii,:,:)), squeeze (ku_face5(3,ii,:,:))); 
     plot3 (squeeze (ku_face6(1,ii,:,:)), squeeze (ku_face6(2,ii,:,:)), squeeze (ku_face6(3,ii,:,:))); 
   end
   for ii = 1:numel(knt2)
     plot3 (squeeze (kv_face1(1,:,ii,:)), squeeze (kv_face1(2,:,ii,:)), squeeze (kv_face1(3,:,ii,:))); 
     plot3 (squeeze (kv_face2(1,:,ii,:)), squeeze (kv_face2(2,:,ii,:)), squeeze (kv_face2(3,:,ii,:))); 
     plot3 (squeeze (kv_face5(1,:,ii,:)), squeeze (kv_face5(2,:,ii,:)), squeeze (kv_face5(3,:,ii,:))); 
     plot3 (squeeze (kv_face6(1,:,ii,:)), squeeze (kv_face6(2,:,ii,:)), squeeze (kv_face6(3,:,ii,:))); 
   end
   for ii = 1:numel(knt3)
     plot3 (squeeze (kw_face1(1,:,:,ii)), squeeze(kw_face1(2,:,:,ii)), squeeze (kw_face1(3,:,:,ii))); 
     plot3 (squeeze (kw_face2(1,:,:,ii)), squeeze(kw_face2(2,:,:,ii)), squeeze (kw_face2(3,:,:,ii))); 
     plot3 (squeeze (kw_face3(1,:,:,ii)), squeeze(kw_face3(2,:,:,ii)), squeeze (kw_face3(3,:,:,ii))); 
     plot3 (squeeze (kw_face4(1,:,:,ii)), squeeze(kw_face4(2,:,:,ii)), squeeze (kw_face4(3,:,:,ii))); 
   end
 end

else % plot a NURBS curve
  nsub = 1000;
  nrbplot (nurbs, nsub);
  hold on

  % And plot the knots
   p = nrbeval (nurbs, unique (nurbs.knots));

   if (any (nurbs.coefs(3,:))) % plot a 3D curve
     plot3 (p(1,:), p(2,:), p(3,:), 'x'); 
   else                     % plot a 2D curve
     plot (p(1,:), p(2,:), 'x'); 
   end

end

if (~hold_flag)
  hold off
end

end
%!demo
%! crv = nrbtestcrv;
%! nrbkntplot(crv)
%! title('Test curve')
%! hold off

%!demo
%! sphere = nrbrevolve(nrbcirc(1,[],0.0,pi),[0.0 0.0 0.0],[1.0 0.0 0.0]);
%! nrbkntplot(sphere);
%! title('Ball and torus - surface construction by revolution');
%! hold on;
%! torus = nrbrevolve(nrbcirc(0.2,[0.9 1.0]),[0.0 0.0 0.0],[1.0 0.0 0.0]);
%! nrbkntplot(torus);
%! hold off

%!demo
%! knots = {[0 0 0 1/2 1 1 1] [0 0 0 1 1 1]...
%!          [0 0 0 1/6 2/6 1/2 1/2 4/6 5/6 1 1 1]};
%!
%! coefs = [-1.0000   -0.9734   -0.7071    1.4290    1.0000    3.4172
%!          0    2.4172         0    0.0148   -2.0000   -1.9734
%!          0    2.0000    4.9623    9.4508    4.0000    2.0000
%!     1.0000    1.0000    0.7071    1.0000    1.0000    1.0000
%!    -0.8536         0   -0.6036    1.9571    1.2071    3.5000
%!     0.3536    2.5000    0.2500    0.5429   -1.7071   -1.0000
%!          0    2.0000    4.4900    8.5444    3.4142    2.0000
%!     0.8536    1.0000    0.6036    1.0000    0.8536    1.0000
%!    -0.3536   -4.0000   -0.2500   -1.2929    1.7071    1.0000
%!     0.8536         0    0.6036   -2.7071   -1.2071   -5.0000
%!          0    2.0000    4.4900   10.0711    3.4142    2.0000
%!     0.8536    1.0000    0.6036    1.0000    0.8536    1.0000
%!          0   -4.0000         0    0.7071    2.0000    5.0000
%!     1.0000    4.0000    0.7071   -0.7071   -1.0000   -5.0000
%!          0    2.0000    4.9623   14.4142    4.0000    2.0000
%!     1.0000    1.0000    0.7071    1.0000    1.0000    1.0000
%!    -2.5000   -4.0000   -1.7678    0.7071    1.0000    5.0000
%!          0    4.0000         0   -0.7071   -3.5000   -5.0000
%!          0    2.0000    6.0418   14.4142    4.0000    2.0000
%!     1.0000    1.0000    0.7071    1.0000    1.0000    1.0000
%!    -2.4379         0   -1.7238    2.7071    1.9527    5.0000
%!     0.9527    4.0000    0.6737    1.2929   -3.4379   -1.0000
%!          0    2.0000    6.6827   10.0711    4.0000    2.0000
%!     1.0000    1.0000    0.7071    1.0000    1.0000    1.0000
%!    -0.9734   -1.0000   -0.6883    0.7071    3.4172    1.0000
%!     2.4172         0    1.7092   -1.4142   -1.9734   -2.0000
%!          0    4.0000    6.6827    4.9623    4.0000         0
%!     1.0000    1.0000    0.7071    0.7071    1.0000    1.0000
%!          0   -0.8536         0    0.8536    3.5000    1.2071
%!     2.5000    0.3536    1.7678   -1.2071   -1.0000   -1.7071
%!          0    3.4142    6.0418    4.4900    4.0000         0
%!     1.0000    0.8536    0.7071    0.6036    1.0000    0.8536
%!    -4.0000   -0.3536   -2.8284    1.2071    1.0000    1.7071
%!          0    0.8536         0   -0.8536   -5.0000   -1.2071
%!          0    3.4142    7.1213    4.4900    4.0000         0
%!     1.0000    0.8536    0.7071    0.6036    1.0000    0.8536
%!    -4.0000         0   -2.8284    1.4142    5.0000    2.0000
%!     4.0000    1.0000    2.8284   -0.7071   -5.0000   -1.0000
%!          0    4.0000   10.1924    4.9623    4.0000         0
%!     1.0000    1.0000    0.7071    0.7071    1.0000    1.0000
%!    -4.0000   -2.5000   -2.8284    0.7071    5.0000    1.0000
%!     4.0000         0    2.8284   -2.4749   -5.0000   -3.5000
%!          0    4.0000   10.1924    6.0418    4.0000         0
%!     1.0000    1.0000    0.7071    0.7071    1.0000    1.0000
%!          0   -2.4379         0    1.3808    5.0000    1.9527
%!     4.0000    0.9527    2.8284   -2.4309   -1.0000   -3.4379
%!          0    4.0000    7.1213    6.6827    4.0000         0
%!     1.0000    1.0000    0.7071    0.7071    1.0000    1.0000
%!    -1.0000   -0.9734    0.2071    2.4163    1.0000    3.4172
%!          0    2.4172   -1.2071   -1.3954   -2.0000   -1.9734
%!     2.0000    4.0000    7.0178    6.6827    2.0000         0
%!     1.0000    1.0000    1.0000    0.7071    1.0000    1.0000
%!    -0.8536         0    0.3536    2.4749    1.2071    3.5000
%!     0.3536    2.5000   -0.8536   -0.7071   -1.7071   -1.0000
%!     1.7071    4.0000    6.3498    6.0418    1.7071         0
%!     0.8536    1.0000    0.8536    0.7071    0.8536    1.0000
%!    -0.3536   -4.0000    0.8536    0.7071    1.7071    1.0000
%!     0.8536         0   -0.3536   -3.5355   -1.2071   -5.0000
%!     1.7071    4.0000    6.3498    7.1213    1.7071         0
%!     0.8536    1.0000    0.8536    0.7071    0.8536    1.0000
%!          0   -4.0000    1.2071    3.5355    2.0000    5.0000
%!     1.0000    4.0000   -0.2071   -3.5355   -1.0000   -5.0000
%!     2.0000    4.0000    7.0178   10.1924    2.0000         0
%!     1.0000    1.0000    1.0000    0.7071    1.0000    1.0000
%!    -2.5000   -4.0000   -0.5429    3.5355    1.0000    5.0000
%!          0    4.0000   -1.9571   -3.5355   -3.5000   -5.0000
%!     2.0000    4.0000    8.5444   10.1924    2.0000         0
%!     1.0000    1.0000    1.0000    0.7071    1.0000    1.0000
%!    -2.4379         0   -0.0355    3.5355    1.9527    5.0000
%!     0.9527    4.0000   -1.4497   -0.7071   -3.4379   -1.0000
%!     2.0000    4.0000    9.4508    7.1213    2.0000         0
%!     1.0000    1.0000    1.0000    0.7071    1.0000    1.0000];
%! coefs = reshape (coefs, 4, 4, 3, 9);
%! horseshoe = nrbmak (coefs, knots);
%! nrbkntplot (horseshoe);
