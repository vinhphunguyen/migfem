function [p,w] = nrbeval(nurbs,tt)
% 
% NRBEVAL: Evaluate a NURBS at parametric points.
% 
% Calling Sequences:
% 
%   [p,w] = nrbeval(crv,ut)
%   [p,w] = nrbeval(srf,{ut,vt})
%   [p,w] = nrbeval(vol,{ut,vt,wt})
%   [p,w] = nrbeval(srf,pts)
% 
% INPUT:
% 
%   crv		: NURBS curve, see nrbmak.
% 
%   srf		: NURBS surface, see nrbmak.
%
%   vol		: NURBS volume, see nrbmak.
% 
%   ut		: Parametric evaluation points along U direction.
%
%   vt		: Parametric evaluation points along V direction.
% 
%   wt		: Parametric evaluation points along W direction.
%
%   pts     : Array of scattered points in parametric domain
% 
% OUTPUT:
%
%   p		: Evaluated points on the NURBS curve, surface or volume as 
% 		Cartesian coordinates (x,y,z). If w is included on the lhs argument
% 		list the points are returned as homogeneous coordinates (wx,wy,wz).
% 
%   w		: Weights of the homogeneous coordinates of the evaluated
% 		points. Note inclusion of this argument changes the type 
% 		of coordinates returned in p (see above).
% 
% Description:
% 
%   Evaluation of NURBS curves, surfaces or volume at parametric points along  
%   the U, V and W directions. Either homogeneous coordinates are returned
%   if the weights are requested in the lhs arguments, or as Cartesian coordinates.
%   This function utilises the 'C' interface bspeval.
% 
% Examples:
% 
%   Evaluate the NURBS circle at twenty points from 0.0 to 1.0
% 
%   nrb = nrbcirc;
%   ut = linspace(0.0,1.0,20);
%   p = nrbeval(nrb,ut);
% 
% See also:
%  
%     bspeval
%
% Copyright (C) 2000 Mark Spink 
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
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

if (nargin < 2)
  error('Not enough input arguments');
end

foption = 1;    % output format 3D cartesian coordinates
if (nargout == 2)
  foption = 0;  % output format 4D homogenous coordinates 
end
   
if (~isstruct(nurbs))
  error('NURBS representation is not structure!');
end

if (~strcmp(nurbs.form,'B-NURBS'))
  error('Not a recognised NURBS representation');
end

if (iscell(nurbs.knots))
  if (size(nurbs.knots,2) == 3)
    %% NURBS structure represents a volume

    num1 = nurbs.number(1);
    num2 = nurbs.number(2);
    num3 = nurbs.number(3);
    degree = nurbs.order-1;

    if (iscell(tt))
      nt1 = numel (tt{1});
      nt2 = numel (tt{2});
      nt3 = numel (tt{3});

      %% evaluate along the w direction
      val = reshape (nurbs.coefs, 4*num1*num2, num3);
      val = bspeval (degree(3), val, nurbs.knots{3}, tt{3});
      val = reshape (val, [4 num1 num2 nt3]);

      %% Evaluate along the v direction
      val = permute (val, [1 2 4 3]);
      val = reshape (val, 4*num1*nt3, num2);
      val = bspeval (degree(2), val, nurbs.knots{2}, tt{2});
      val = reshape (val, [4 num1 nt3 nt2]);
      val = permute (val, [1 2 4 3]);

      %% Evaluate along the u direction
      val = permute (val, [1 3 4 2]);
      val = reshape (val, 4*nt2*nt3, num1);
      val = bspeval (degree(1), val, nurbs.knots{1}, tt{1});
      val = reshape (val, [4 nt2 nt3 nt1]);
      val = permute (val, [1 4 2 3]);
      pnts = val;

      p = pnts(1:3,:,:,:);
      w = pnts(4,:,:,:);
      if (foption)
        p = p./repmat(w,[3 1 1 1]);
      end

    else

      %% Evaluate at scattered points
      %% tt(1,:) represents the u direction
      %% tt(2,:) represents the v direction
      %% tt(3,:) represents the w direction

      %% evaluate along the w direction
      nt = size(tt,2);
      val = reshape(nurbs.coefs,4*num1*num2,num3);
      val = bspeval(degree(3),val,nurbs.knots{3},tt(3,:));
      val = reshape(val,[4 num1 num2 nt]);

      %% evaluate along the v direction
      val2 = zeros(4*num1,nt);
      for v = 1:nt
        coefs = reshape(val(:,:,:,v),4*num1,num2);
        val2(:,v) = bspeval(degree(2),coefs,nurbs.knots{2},tt(2,v));
      end
      val2 = reshape(val2,[4 num1 nt]);

      %% evaluate along the u direction
      pnts = zeros(4,nt);
      for v = 1:nt
        coefs = reshape (val2(:,:,v), [4 num1]);
        pnts(:,v) = bspeval(degree(1),coefs,nurbs.knots{1},tt(1,v));
      end

      w = pnts(4,:);
      p = pnts(1:3,:);
      if (foption)
        p = p./repmat(w,[3, 1]);
      end
    end

  elseif (size(nurbs.knots,2) == 2)
    %% NURBS structure represents a surface
  
    num1 = nurbs.number(1);
    num2 = nurbs.number(2);
    degree = nurbs.order-1;

    if (iscell(tt))
      %% Evaluate over a [u,v] grid
      %% tt{1} represents the u direction
      %% tt{2} represents the v direction

      nt1 = length(tt{1});
      nt2 = length(tt{2});
    
      %% Evaluate along the v direction
      val = reshape(nurbs.coefs,4*num1,num2);
      val = bspeval(degree(2),val,nurbs.knots{2},tt{2});
      val = reshape(val,[4 num1 nt2]);
    
      %% Evaluate along the u direction
      val = permute(val,[1 3 2]);
      val = reshape(val,4*nt2,num1);
      val = bspeval(degree(1),val,nurbs.knots{1},tt{1});
      val = reshape(val,[4 nt2 nt1]);
      val = permute(val,[1 3 2]);

      w = val(4,:,:);
      p = val(1:3,:,:);
      if (foption)
	p = p./repmat(w,[3 1 1]);
      end

    else

      %% Evaluate at scattered points
      %% tt(1,:) represents the u direction
      %% tt(2,:) represents the v direction

      nt = size(tt,2);

      val = reshape(nurbs.coefs,4*num1,num2);
      val = bspeval(degree(2),val,nurbs.knots{2},tt(2,:));
      val = reshape(val,[4 num1 nt]);


      %% evaluate along the u direction
      pnts = zeros(4,nt);
      for v = 1:nt
	coefs = reshape (val(:,:,v), [4 num1]);
	pnts(:,v) = bspeval(degree(1),coefs,nurbs.knots{1},tt(1,v));
      end

      w = pnts(4,:);
      p = pnts(1:3,:);
      if (foption)
	p = p./repmat(w,[3, 1]);
      end
        
    end

  end
else

  %% NURBS structure represents a curve
  %%  tt represent a vector of parametric points in the u direction
  
  val = bspeval(nurbs.order-1,nurbs.coefs,nurbs.knots,tt);   

  w = val(4,:);
  p = val(1:3,:);
  if foption
    p = p./repmat(w,3,1);
  end

end

end

%!demo
%! srf = nrbtestsrf;
%! p = nrbeval(srf,{linspace(0.0,1.0,20) linspace(0.0,1.0,20)});
%! h = surf(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)));
%! title('Test surface.');
%! hold off

%!test
%! knots{1} = [0 0 0 1 1 1];
%! knots{2} = [0 0 0 .5 1 1 1];
%! knots{3} = [0 0 0 0 1 1 1 1];
%! cx = [0 0.5 1]; nx = length(cx);
%! cy = [0 0.25 0.75 1]; ny = length(cy);
%! cz = [0 1/3 2/3 1]; nz = length(cz);
%! coefs(1,:,:,:) = repmat(reshape(cx,nx,1,1),[1 ny nz]);
%! coefs(2,:,:,:) = repmat(reshape(cy,1,ny,1),[nx 1 nz]);
%! coefs(3,:,:,:) = repmat(reshape(cz,1,1,nz),[nx ny 1]);
%! coefs(4,:,:,:) = 1;
%! nurbs = nrbmak(coefs, knots);
%! x = rand(5,1); y = rand(5,1); z = rand(5,1);
%! tt = [x y z]';
%! points = nrbeval(nurbs,tt);
%!
%! assert(points,tt,1e-10)
%!
%!test
%! knots{1} = [0 0 0 1 1 1];
%! knots{2} = [0 0 0 0 1 1 1 1];
%! knots{3} = [0 0 1 1];
%! cx = [0 0 1]; nx = length(cx);
%! cy = [0 0 0 1]; ny = length(cy);
%! cz = [0 1]; nz = length(cz);
%! coefs(1,:,:,:) = repmat(reshape(cx,nx,1,1),[1 ny nz]);
%! coefs(2,:,:,:) = repmat(reshape(cy,1,ny,1),[nx 1 nz]);
%! coefs(3,:,:,:) = repmat(reshape(cz,1,1,nz),[nx ny 1]);
%! coefs(4,:,:,:) = 1;
%! nurbs = nrbmak(coefs, knots);
%! x = rand(5,1); y = rand(5,1); z = rand(5,1);
%! tt = [x y z]';
%! points = nrbeval(nurbs,tt);
%! assert(points,[x.^2 y.^3 z]',1e-10);
%!
%!test
%! knots{1} = [0 0 0 1 1 1];
%! knots{2} = [0 0 0 0 1 1 1 1];
%! knots{3} = [0 0 1 1];
%! cx = [0 0 1]; nx = length(cx);
%! cy = [0 0 0 1]; ny = length(cy);
%! cz = [0 1]; nz = length(cz);
%! coefs(1,:,:,:) = repmat(reshape(cx,nx,1,1),[1 ny nz]);
%! coefs(2,:,:,:) = repmat(reshape(cy,1,ny,1),[nx 1 nz]);
%! coefs(3,:,:,:) = repmat(reshape(cz,1,1,nz),[nx ny 1]);
%! coefs(4,:,:,:) = 1;
%! coefs = coefs([2 1 3 4],:,:,:);
%! nurbs = nrbmak(coefs, knots);
%! x = rand(5,1); y = rand(5,1); z = rand(5,1);
%! tt = [x y z]';
%! points = nrbeval(nurbs,tt);
%! [y.^3 x.^2 z]';
%! assert(points,[y.^3 x.^2 z]',1e-10);
