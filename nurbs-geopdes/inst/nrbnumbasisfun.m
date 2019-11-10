function idx = nrbnumbasisfun (points, nrb)
%
% NRBNUMBASISFUN:  Numbering of basis functions for NURBS
%
% Calling Sequence:
% 
%   N      = nrbnumbasisfun (u, crv)
%   N      = nrbnumbasisfun ({u, v}, srf)
%   N      = nrbnumbasisfun (p, srf)
%
%    INPUT:
%   
%      u or p(1,:,:)  - parametric points along u direction
%      v or p(2,:,:)  - parametric points along v direction
%      crv - NURBS curve
%      srf - NURBS surface
%   
%    OUTPUT:
%
%      N - Indices of the basis functions that are nonvanishing at each
%          point. size(N) == size(B)
%
%    Copyright (C) 2009 Carlo de Falco
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

  if (   (nargin<2) ...
      || (nargout>1) ...
      || (~isstruct(nrb)) ...
      || (iscell(points) && ~iscell(nrb.knots)) ...
      || (~iscell(points) && iscell(nrb.knots) && (size(points,1)~=2)) ...
      )
    error('Incorrect input arguments in nrbnumbasisfun');
  end


  if (~iscell(nrb.knots))          %% NURBS curve
    
    iv  = findspan (nrb.number-1, nrb.order-1, points, nrb.knots);
    idx = numbasisfun (iv, points, nrb.order-1, nrb.knots);
    
  elseif size(nrb.knots,2) == 2  %% NURBS surface

    if (iscell(points))
      [v, u] = meshgrid(points{2}, points{1});
      p(1,:,:) = u;
      p(2,:,:) = v;
      p = reshape(p, 2, []);
    else
      p = points;
    end
    
    idx = nrb_srf_numbasisfun__ (p, nrb); 
  else
    error('The function nrbnumbasisfun is not yet ready for volumes')      
  end
  
end


%!test
%! p = 2;   q = 3;   m = 4; n = 5;
%! Lx  = 1; Ly  = 1; 
%! nrb = nrb4surf   ([0 0], [1 0], [0 1], [1 1]);
%! nrb = nrbdegelev (nrb, [p-1, q-1]);
%! ikx = linspace(0,1,m); iky = linspace(0,1,n);
%! nrb = nrbkntins  (nrb, {ikx(2:end-1), iky(2:end-1)});
%! nrb.coefs (4,:,:) = nrb.coefs (4,:,:) + rand (size (nrb.coefs (4,:,:)));
%! u = rand (1, 30); v = rand (1, 10);
%! u = (u-min (u))/max (u-min (u));
%! v = (v-min (v))/max (v-min (v));
%! N = nrbnumbasisfun ({u, v}, nrb);
%! assert (all (all (N>0)), true)
%! assert (all (all (N <= prod (nrb.number))), true)
%! assert (max (max (N)), prod (nrb.number))
%! assert (min (min (N)), 1)