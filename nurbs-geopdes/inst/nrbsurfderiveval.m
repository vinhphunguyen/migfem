function skl = nrbsurfderiveval (srf, uv, d) 

%
% NRBSURFDERIVEVAL: Evaluate n-th order derivatives of a NURBS surface
%
% usage: skl = nrbsurfderiveval (srf, [u; v], d) 
%
%   INPUT:
%
%   srf   : NURBS surface structure, see nrbmak
%
%   u, v  : parametric coordinates of the point where we compute the
%      derivatives
%
%   d     : number of partial derivatives to compute
%
%
%   OUTPUT: 
%
%   skl (i, j, k, l) = i-th component derived j-1,k-1 times at the
%     l-th point.
%
% Adaptation of algorithm A4.4 from the NURBS book, pg137
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
  
 skl = zeros (3, d+1, d+1, size (uv, 2));
 
 for iu = 1:size(uv, 2);
    wders = squeeze  (surfderiveval (srf.number(1)-1, srf.order(1)-1,  ...
				     srf.knots{1}, srf.number(2)-1,  ...
				     srf.order(2)-1, srf.knots{2},  ...
				     squeeze (srf.coefs(4, :, :)), uv(1,iu), ...
				     uv(2,iu), d));      

  for idim = 1:3
      Aders = squeeze  (surfderiveval (srf.number(1)-1, srf.order(1)-1,  ...
			srf.knots{1}, srf.number(2)-1,  ...
			srf.order(2)-1, srf.knots{2},  ...
			squeeze (srf.coefs(idim, :, :)), uv(1,iu),...
			uv(2,iu), d));      
      
   for k=0:d
    for l=0:d-k
	  v = Aders(k+1, l+1);
      for j=1:l
	    v = v - nchoosek(l,j)*wders(1,j+1)*skl(idim, k+1, l-j+1,iu);
      end
      for i=1:k
	    v = v - nchoosek(k,i)*wders(i+1,1)*skl(idim, k-i+1, l+1, iu);
	    v2 =0;
        for j=1:l
	      v2 = v2 + nchoosek(l,j)*wders(i+1,j+1)*skl(idim, k-i+1, l-j+1, iu);
        end
	    v = v - nchoosek(k,i)*v2;
      end
	  skl(idim, k+1, l+1, iu) = v/wders(1,1);
    end
   end
  end
    
 end

end

%!test
%! k = [0 0 1 1];
%! c = [0 1];
%! [coef(2,:,:), coef(1,:,:)] = meshgrid (c, c);
%! coef(3,:,:) = coef(1,:,:);
%! srf = nrbmak (coef, {k, k});
%! [u, v] = meshgrid (linspace(0,1,11));
%! uv = [u(:)';v(:)'];
%! skl = nrbsurfderiveval (srf, uv, 0);
%! aux = nrbeval(srf,uv);
%! assert (squeeze (skl (1:2,1,1,:)), aux(1:2,:), 1e3*eps)


%!test
%! k = [0 0 1 1];
%! c = [0 1];
%! [coef(2,:,:), coef(1,:,:)] = meshgrid (c, c);
%! coef(3,:,:) = coef(1,:,:);
%! srf = nrbmak (coef, {k, k});
%! srf = nrbkntins (srf, {[], rand(2,1)});
%! [u, v] = meshgrid (linspace(0,1,11));
%! uv = [u(:)';v(:)'];
%! skl = nrbsurfderiveval (srf, uv, 0);
%! aux = nrbeval(srf,uv);
%! assert (squeeze (skl (1:2,1,1,:)), aux(1:2,:), 1e3*eps)

%!shared srf, uv
%!test 
%! k = [0 0 0 1 1 1];
%! c = [0 1/2 1];
%! [coef(1,:,:), coef(2,:,:)] = meshgrid (c, c);
%! coef(3,:,:) = coef(1,:,:);
%! srf = nrbmak (coef, {k, k});
%! ders= nrbderiv (srf);
%! [u, v] = meshgrid (linspace(0,1,11));
%! uv = [u(:)';v(:)'];
%! skl = nrbsurfderiveval (srf, uv, 1);
%! [fun, der] = nrbdeval (srf, ders, uv);
%! assert (squeeze (skl (1:2,1,1,:)), fun(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,2,1,:)), der{1}(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,1,2,:)), der{2}(1:2,:), 1e3*eps)
%!
%!test 
%! srf = nrbdegelev (srf, [3, 1]);
%! ders= nrbderiv (srf);
%! [fun, der] = nrbdeval (srf, ders, uv);
%! skl = nrbsurfderiveval (srf, uv, 1);
%! assert (squeeze (skl (1:2,1,1,:)), fun(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,2,1,:)), der{1}(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,1,2,:)), der{2}(1:2,:), 1e3*eps)

%!shared uv
%!test 
%! k = [0 0 0 1 1 1];
%! c = [0 1/2 1];
%! [coef(2,:,:), coef(1,:,:)] = meshgrid (c, c);
%! coef(3,:,:) = coef(1,:,:);
%! srf = nrbmak (coef, {k, k});
%! ders= nrbderiv (srf);
%! [u, v] = meshgrid (linspace(0,1,11));
%! uv = [u(:)';v(:)'];
%! skl = nrbsurfderiveval (srf, uv, 1);
%! [fun, der] = nrbdeval (srf, ders, uv);
%! assert (squeeze (skl (1:2,1,1,:)), fun(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,2,1,:)), der{1}(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,1,2,:)), der{2}(1:2,:), 1e3*eps)
%!
%!test 
%! p = 3;  q = 3;
%! mcp = 5; ncp = 5;
%! Lx  = 10*rand(1); Ly  = Lx;
%! srf = nrbdegelev (nrb4surf ([0 0], [Lx, 0], [0 Ly], [Lx Ly]), [p-1, q-1]);
%! %%srf = nrbkntins (srf, {linspace(0,1,mcp-p+2)(2:end-1), linspace(0,1,ncp-q+2)(2:end-1)});
%! %%srf.coefs = permute (srf.coefs, [1 3 2]);
%! ders= nrbderiv (srf);
%! [fun, der] = nrbdeval (srf, ders, uv);
%! skl = nrbsurfderiveval (srf, uv, 1);
%! assert (squeeze (skl (1:2,1,1,:)), fun(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,2,1,:)), der{1}(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,1,2,:)), der{2}(1:2,:), 1e3*eps)

%!shared srf, uv, P, dPdx, d2Pdx2, c1, c2
%!test
%! [u, v] = meshgrid (linspace(0,1,10));
%! uv = [u(:)';v(:)'];
%! c1 = nrbmak([0 1/2 1; 0 1 0],[0 0 0 1 1 1]);
%! c1 = nrbtform (c1, vecrotx (pi/2));
%! c2  = nrbtform(c1, vectrans([0 1 0]));
%! srf = nrbdegelev (nrbruled (c1, c2), [3, 1]);
%! skl = nrbsurfderiveval (srf, uv, 2);
%! P = squeeze(skl(:,1,1,:));
%! dPdx = squeeze(skl(:,2,1,:));
%! d2Pdx2 = squeeze(skl(:,3,1,:));
%!assert(P(3,:), 2*(P(1,:)-P(1,:).^2),100*eps)
%!assert(dPdx(3,:), 2-4*P(1,:), 100*eps)
%!assert(d2Pdx2(3,:), -4+0*P(1,:), 100*eps)
%! srf = nrbdegelev (nrbruled (c1, c2), [5, 6]);
%! skl = nrbsurfderiveval (srf, uv, 2);
%! P = squeeze(skl(:,1,1,:));
%! dPdx = squeeze(skl(:,2,1,:));
%! d2Pdx2 = squeeze(skl(:,3,1,:));
%! aux = nrbeval(srf,uv);
%! assert (squeeze (skl (1:2,1,1,:)), aux(1:2,:), 1e3*eps)
%!assert(P(3,:), 2*(P(1,:)-P(1,:).^2),100*eps)
%!assert(dPdx(3,:), 2-4*P(1,:), 100*eps)
%!assert(d2Pdx2(3,:), -4+0*P(1,:), 100*eps)
%!
%!test
%! skl = nrbsurfderiveval (srf, uv, 0);
%! aux = nrbeval (srf, uv);
%! assert (squeeze (skl (1:2,1,1,:)), aux(1:2,:), 1e3*eps)

%!shared dPdu, d2Pdu2, P, srf, uv
%!test
%! [u, v] = meshgrid (linspace(0,1,10));
%! uv = [u(:)';v(:)'];
%! c1 = nrbmak([0 1/2 1; 0.1 1.6 1.1; 0 0 0],[0 0 0 1 1 1]);
%! c2 = nrbmak([0 1/2 1; 0.1 1.6 1.1; 1 1 1],[0 0 0 1 1 1]);
%! srf = nrbdegelev (nrbruled (c1, c2), [0, 1]);
%! skl = nrbsurfderiveval (srf, uv, 2);
%! P = squeeze(skl(:,1,1,:));
%! dPdu = squeeze(skl(:,2,1,:));
%! dPdv = squeeze(skl(:,1,2,:));
%! d2Pdu2 = squeeze(skl(:,3,1,:));
%! aux = nrbeval(srf,uv);
%! assert (squeeze (skl (1:2,1,1,:)), aux(1:2,:), 1e3*eps)
%!assert(dPdu(2,:), 3-4*P(1,:),100*eps)
%!assert(d2Pdu2(2,:), -4+0*P(1,:),100*eps)
%!
%!test
%! skl = nrbsurfderiveval (srf, uv, 0);
%! aux = nrbeval(srf,uv);
%! assert (squeeze (skl (1:2,1,1,:)), aux(1:2,:), 1e3*eps)

%!test
%! srf = nrb4surf([0 0], [1 0], [0 1], [1 1]);
%! geo = nrbdegelev (srf, [3 3]);
%1 geo.coefs (4, 2:end-1, 2:end-1) = geo.coefs (4, 2:end-1, 2:end-1) + .1 * rand (1, geo.number(1)-2, geo.number(2)-2);
%! geo = nrbkntins (geo, {[.1:.1:.9], [.2:.2:.8]});
%! [u, v] = meshgrid (linspace(0,1,10));
%! uv = [u(:)';v(:)'];
%! skl = nrbsurfderiveval (geo, uv, 2);
%! dgeo = nrbderiv (geo);
%! [pnts, ders] = nrbdeval (geo, dgeo, uv);
%! assert (ders{1}, squeeze(skl(:,2,1,:)), 1e-9)
%! assert (ders{2}, squeeze(skl(:,1,2,:)), 1e-9)

%!test
%! crv = nrbline ([1 0], [2 0]);
%! srf = nrbrevolve (crv, [0 0 0], [0 0 1], pi/2);
%! srf = nrbtransp (srf);
%! [v, u] = meshgrid (linspace (0, 1, 11));
%! uv = [u(:)'; v(:)'];
%! skl = nrbsurfderiveval (srf, uv, 2);
%! c = sqrt(2);
%! w      = @(x, y) (2 - c)*y.^2 + (c-2)*y + 1;
%! dwdy   = @(x, y) 2*(2-c)*y + c - 2;
%! d2wdy2 = @(x, y) 2*(2-c);
%! F1 = @(x, y) (x+1) .* ((1-y).^2 + c*y.*(1-y)) ./ w(x,y);
%! F2 = @(x, y) (x+1) .* (y.^2 + c*y.*(1-y)) ./ w(x,y);
%! dF1dx = @(x, y) ((1-y).^2 + c*y.*(1-y)) ./ w(x,y);
%! dF2dx = @(x, y) (y.^2 + c*y.*(1-y)) ./ w(x,y);
%! dF1dy = @(x, y) (x+1) .* ((2 - 2*c)*y + c - 2) ./ w(x,y) - (x+1) .* ((1-y).^2 + c*y.*(1-y)) .* dwdy(x,y) ./ w(x,y).^2;
%! dF2dy = @(x, y) (x+1) .* ((2 - 2*c)*y + c) ./ w(x,y) - (x+1) .* (y.^2 + c*y.*(1-y)) .* dwdy(x,y) ./ w(x,y).^2;
%! d2F1dx2 = @(x, y) zeros (size (x));
%! d2F2dx2 = @(x, y) zeros (size (x));
%! d2F1dxdy = @(x, y) ((2 - 2*c)*y + c - 2) ./ w(x,y) - ((1-y).^2 + c*y.*(1-y)) .* dwdy(x,y) ./ w(x,y).^2;
%! d2F2dxdy = @(x, y) ((2 - 2*c)*y + c) ./ w(x,y) - (y.^2 + c*y.*(1-y)) .* dwdy(x,y) ./ w(x,y).^2;
%! d2F1dy2  = @(x, y) (x+1)*(2 - 2*c) ./ w(x,y) - 2*(x+1) .* ((2 - 2*c)*y + c - 2) .* dwdy(x,y) ./ w(x,y).^2 - ...
%!                    (x+1) .* ((1-y).^2 + c*y.*(1-y)) * d2wdy2(x,y) ./ w(x,y).^2 + ...
%!                    2 * (x+1) .* ((1-y).^2 + c*y.*(1-y)) .* w(x,y) .*dwdy(x,y).^2 ./ w(x,y).^4;
%! d2F2dy2  = @(x, y) (x+1)*(2 - 2*c) ./ w(x,y) - 2*(x+1) .* ((2 - 2*c)*y + c) .* dwdy(x,y) ./ w(x,y).^2 - ...
%!                    (x+1) .* (y.^2 + c*y.*(1-y)) * d2wdy2(x,y) ./ w(x,y).^2 + ...
%!                    2 * (x+1) .* (y.^2 + c*y.*(1-y)) .* w(x,y) .*dwdy(x,y).^2 ./ w(x,y).^4;
%! assert ([F1(u(:),v(:)), F2(u(:),v(:))], squeeze(skl(1:2,1,1,:))', 1e2*eps);
%! assert ([dF1dx(u(:),v(:)), dF2dx(u(:),v(:))], squeeze(skl(1:2,2,1,:))', 1e2*eps);
%! assert ([dF1dy(u(:),v(:)), dF2dy(u(:),v(:))], squeeze(skl(1:2,1,2,:))', 1e2*eps);
%! assert ([d2F1dx2(u(:),v(:)), d2F2dx2(u(:),v(:))], squeeze(skl(1:2,3,1,:))', 1e2*eps);
%! assert ([d2F1dxdy(u(:),v(:)), d2F2dxdy(u(:),v(:))], squeeze(skl(1:2,2,2,:))', 1e2*eps);
%! assert ([d2F1dy2(u(:),v(:)), d2F2dy2(u(:),v(:))], squeeze(skl(1:2,1,3,:))', 1e2*eps);

%!test
%! knots = {[0 0 1 1] [0 0 1 1]};
%! coefs(:,1,1) = [0;0;0;1];
%! coefs(:,2,1) = [1;0;0;1];
%! coefs(:,1,2) = [0;1;0;1];
%! coefs(:,2,2) = [1;1;1;2];
%! srf = nrbmak (coefs, knots);
%! [v, u] = meshgrid (linspace (0, 1, 3));
%! uv = [u(:)'; v(:)'];
%! skl = nrbsurfderiveval (srf, uv, 2);
%! w = @(x, y) x.*y + 1;
%! F1 = @(x, y) x ./ w(x,y);
%! F2 = @(x, y) y ./ w(x,y);
%! F3 = @(x, y) x .* y ./ w(x,y);
%! dF1dx = @(x, y) 1./w(x,y) - x.*y./w(x,y).^2;
%! dF1dy = @(x, y)  - x.^2./w(x,y).^2;
%! dF2dx = @(x, y)  - y.^2./w(x,y).^2;
%! dF2dy = @(x, y) 1./w(x,y) - x.*y./w(x,y).^2;
%! dF3dx = @(x, y) y./w(x,y) - x.*(y./w(x,y)).^2;
%! dF3dy = @(x, y) x./w(x,y) - y.*(x./w(x,y)).^2;
%! d2F1dx2  = @(x, y) -2*y./w(x,y).^2 + 2*x.*y.^2./w(x,y).^3;
%! d2F1dy2  = @(x, y) 2*x.^3./w(x,y).^3;
%! d2F1dxdy = @(x, y) -x./w(x,y).^2 - x./w(x,y).^2 + 2*x.^2.*y./w(x,y).^3;
%! d2F2dx2  = @(x, y) 2*y.^3./w(x,y).^3;
%! d2F2dy2  = @(x, y) -2*x./w(x,y).^2 + 2*y.*x.^2./w(x,y).^3;
%! d2F2dxdy = @(x, y) -y./w(x,y).^2 - y./w(x,y).^2 + 2*y.^2.*x./w(x,y).^3;
%! d2F3dx2  = @(x, y) -2*y.^2./w(x,y).^2 + 2*x.*y.^3./w(x,y).^3;
%! d2F3dy2  = @(x, y) -2*x.^2./w(x,y).^2 + 2*y.*x.^3./w(x,y).^3;
%! d2F3dxdy = @(x, y) 1./w(x,y) - 3*x.*y./w(x,y).^2 + 2*(x.*y).^2./w(x,y).^3;
%! assert ([F1(u(:),v(:)), F2(u(:),v(:)), F3(u(:),v(:))], squeeze(skl(1:3,1,1,:))', 1e2*eps);
%! assert ([dF1dx(u(:),v(:)), dF2dx(u(:),v(:)), dF3dx(u(:),v(:))], squeeze(skl(1:3,2,1,:))', 1e2*eps);
%! assert ([dF1dy(u(:),v(:)), dF2dy(u(:),v(:)), dF3dy(u(:),v(:))], squeeze(skl(1:3,1,2,:))', 1e2*eps);
%! assert ([d2F1dx2(u(:),v(:)), d2F2dx2(u(:),v(:)), d2F3dx2(u(:),v(:))], squeeze(skl(1:3,3,1,:))', 1e2*eps);
%! assert ([d2F1dy2(u(:),v(:)), d2F2dy2(u(:),v(:)), d2F3dy2(u(:),v(:))], squeeze(skl(1:3,1,3,:))', 1e2*eps);
%! assert ([d2F1dxdy(u(:),v(:)), d2F2dxdy(u(:),v(:)), d2F3dxdy(u(:),v(:))], squeeze(skl(1:3,2,2,:))', 1e2*eps);
