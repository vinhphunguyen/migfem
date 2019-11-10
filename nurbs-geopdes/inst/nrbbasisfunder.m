function varargout = nrbbasisfunder (points, nrb)

% NRBBASISFUNDER:  NURBS basis functions derivatives
%
% Calling Sequence:
% 
%   Bu          = nrbbasisfunder (u, crv)
%   [Bu, N]     = nrbbasisfunder (u, crv)
%   [Bu, Bv]    = nrbbasisfunder ({u, v}, srf)
%   [Bu, Bv, N] = nrbbasisfunder ({u, v}, srf)
%   [Bu, Bv, N] = nrbbasisfunder (p, srf)
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
%      Bu - Basis functions derivatives WRT direction u
%           size(Bu)=[numel(u),(p+1)] for curves
%           or [numel(u)*numel(v), (p+1)*(q+1)] for surfaces
%
%      Bv - Basis functions derivatives WRT direction v
%           size(Bv)=[numel(v),(p+1)] for curves
%           or [numel(u)*numel(v), (p+1)*(q+1)] for surfaces
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
      || (nargout>3) ...
      || (~isstruct(nrb)) ...
      || (iscell(points) && ~iscell(nrb.knots)) ...
      || (~iscell(points) && iscell(nrb.knots) && (size(points,1)~=2)) ...
      || (~iscell(nrb.knots) && (nargout>2)) ...
      )
    error('Incorrect input arguments in nrbbasisfun');
  end
                            
  if (~iscell(nrb.knots))         %% NURBS curve
    
    [varargout{1}, varargout{2}] = nrb_crv_basisfun_der__ (points, nrb);

  elseif size(nrb.knots,2) == 2 %% NURBS surface

    if (iscell(points))
      [v, u] = meshgrid(points{2}, points{1});
      p = [u(:), v(:)]';
    else
      p = points;
    end
    
    [varargout{1}, varargout{2}, varargout{3}] = nrb_srf_basisfun_der__ (p, nrb);

  else                            %% NURBS volume
    error('The function nrbbasisfunder is not yet ready for volumes')
  end
end
  
%!demo
%! U = [0 0 0 0 1 1 1 1];
%! x = [0 1/3 2/3 1] ;
%! y = [0 0 0 0];
%! w = [1 1 1 1];
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 30);
%! [Bu, id] = nrbbasisfunder (u, nrb);
%! plot(u, Bu)
%! title('Derivatives of the cubic Bernstein polynomials')
%! hold off

%!test
%! U = [0 0 0 0 1 1 1 1];
%! x = [0 1/3 2/3 1] ;
%! y = [0 0 0 0];
%! w = rand(1,4);
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 30);
%! [Bu, id] = nrbbasisfunder (u, nrb);
%! #plot(u, Bu)
%! assert (sum(Bu, 2), zeros(numel(u), 1), 1e-10), 

%!test
%! U = [0 0 0 0 1/2 1 1 1 1];
%! x = [0 1/4 1/2 3/4 1] ;
%! y = [0 0 0 0 0];
%! w = rand(1,5);
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 300); 
%! [Bu, id] = nrbbasisfunder (u, nrb); 
%! assert (sum(Bu, 2), zeros(numel(u), 1), 1e-10)

%!test
%! p = 2;   q = 3;   m = 4; n = 5;
%! Lx  = 1; Ly  = 1; 
%! nrb = nrb4surf   ([0 0], [1 0], [0 1], [1 1]);
%! nrb = nrbdegelev (nrb, [p-1, q-1]);
%! aux1 = linspace(0,1,m); aux2 = linspace(0,1,n);
%! nrb = nrbkntins  (nrb, {aux1(2:end-1), aux2(2:end-1)});
%! nrb.coefs (4,:,:) = nrb.coefs(4,:,:) + rand (size (nrb.coefs (4,:,:)));
%! [Bu, Bv, N] = nrbbasisfunder ({rand(1, 20), rand(1, 20)}, nrb);
%! #plot3(squeeze(u(1,:,:)), squeeze(u(2,:,:)), reshape(Bu(:,10), 20, 20),'o')
%! assert (sum (Bu, 2), zeros(20^2, 1), 1e-10)

