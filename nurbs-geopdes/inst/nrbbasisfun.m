function [B, id] = nrbbasisfun (points, nrb)

% NRBBASISFUN: Basis functions for NURBS
%
% Calling Sequence:
% 
%    B     = nrbbasisfun (u, crv)
%    B     = nrbbasisfun ({u, v}, srf)
%   [B, N] = nrbbasisfun ({u, v}, srf)
%   [B, N] = nrbbasisfun (p, srf)
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
%      B - Value of the basis functions at the points
%          size(B)=[numel(u),(p+1)] for curves
%          or [numel(u)*numel(v), (p+1)*(q+1)] for surfaces
%
%      N - Indices of the basis functions that are nonvanishing at each
%          point. size(N) == size(B)
%   
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
      || (nargout>2) ...
      || (~isstruct(nrb)) ...
      || (iscell(points) && ~iscell(nrb.knots)) ...
      || (~iscell(points) && iscell(nrb.knots) && (size(points,1)~=2)) ...
      || (~iscell(nrb.knots) && (nargout>1)) ...
      )
    error('Incorrect input arguments in nrbbasisfun');
  end
                            
  if (~iscell(nrb.knots))         %% NURBS curve
    
    [B, id] = nrb_crv_basisfun__ (points, nrb);
    
  elseif size(nrb.knots,2) == 2 %% NURBS surface
    if (iscell(points))
      [v, u] = meshgrid(points{2}, points{1});
      p = [u(:), v(:)]';
    else
      p = points;
    end
    
    [B, id] = nrb_srf_basisfun__ (p, nrb); 

  else                            %% NURBS volume
    error('The function nrbbasisfun is not yet ready for volumes')
  end
end  

%!demo
%! U = [0 0 0 0 1 1 1 1];
%! x = [0 1/3 2/3 1] ;
%! y = [0 0 0 0];
%! w = [1 1 1 1];
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 30);
%! B = nrbbasisfun (u, nrb);
%! xplot = sum(bsxfun(@(x,y) x.*y, B, x),2);
%! plot(xplot, B)
%! title('Cubic Bernstein polynomials')
%! hold off

%!test
%! U = [0 0 0 0 1 1 1 1];
%! x = [0 1/3 2/3 1] ;
%! y = [0 0 0 0];
%! w = rand(1,4);
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 30);
%! B = nrbbasisfun (u, nrb);
%! xplot = sum(bsxfun(@(x,y) x.*y, B, x),2);
%!
%! yy = y; yy(1) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U); 
%! aux = nrbeval(nrb2,u);
%! %figure, plot(xplot, B(:,1), aux(1,:).', w(1)*aux(2,:).')
%! assert(B(:,1), w(1)*aux(2,:).', 1e-6)
%! 
%! yy = y; yy(2) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U);
%! aux = nrbeval(nrb2, u);
%! %figure, plot(xplot, B(:,2), aux(1,:).', w(2)*aux(2,:).')
%! assert(B(:,2), w(2)*aux(2,:).', 1e-6)
%!
%! yy = y; yy(3) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U);
%! aux = nrbeval(nrb2,u);
%! %figure, plot(xplot, B(:,3), aux(1,:).', w(3)*aux(2,:).')
%! assert(B(:,3), w(3)*aux(2,:).', 1e-6)
%!
%! yy = y; yy(4) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U);
%! aux = nrbeval(nrb2,u);
%! %figure, plot(xplot, B(:,4), aux(1,:).', w(4)*aux(2,:).')
%! assert(B(:,4), w(4)*aux(2,:).', 1e-6)

%!test
%! p = 2;   q = 3;   m = 4; n = 5;
%! Lx  = 1; Ly  = 1; 
%! nrb = nrb4surf   ([0 0], [1 0], [0 1], [1 1]);
%! nrb = nrbdegelev (nrb, [p-1, q-1]);
%! aux1 = linspace(0,1,m); aux2 = linspace(0,1,n);
%! nrb = nrbkntins  (nrb, {aux1(2:end-1), aux2(2:end-1)});
%! u = rand (1, 30); v = rand (1, 10);
%! [B, N] = nrbbasisfun ({u, v}, nrb);
%! assert (sum(B, 2), ones(300, 1), 1e-6)
%! assert (all (all (B<=1)), true)
%! assert (all (all (B>=0)), true)
%! assert (all (all (N>0)), true)
%! assert (all (all (N <= prod (nrb.number))), true)
%! assert (max (max (N)),prod (nrb.number))
%! assert (min (min (N)),1)