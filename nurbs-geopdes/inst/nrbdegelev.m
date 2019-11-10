function inurbs = nrbdegelev(nurbs, ntimes)
% 
% NRBDEGELEV: Elevate the degree of the NURBS curve, surface or volume.
% 
% Calling Sequence:
% 
%   ecrv = nrbdegelev(crv,utimes);
%   esrf = nrbdegelev(srf,[utimes,vtimes]);
%   evol = nrbdegelev(vol,[utimes,vtimes,wtimes]);
% 
% INPUT:
% 
%   crv		: NURBS curve, see nrbmak.
% 
%   srf		: NURBS surface, see nrbmak.
% 
%   vol		: NURBS volume, see nrbmak.
% 
%   utimes	: Increase the degree along U direction utimes.
% 
%   vtimes	: Increase the degree along V direction vtimes.
% 
%   wtimes	: Increase the degree along W direction vtimes.
%
% OUTPUT:
%
%   ecrv	: new NURBS structure for a curve with degree elevated.
% 
%   esrf	: new NURBS structure for a surface with degree elevated.
% 
%   evol	: new NURBS structure for a volume with degree elevated.
% 
% 
% Description:
% 
%   Degree elevates the NURBS curve or surface. This function uses the
%   B-Spline function bspdegelev, which interface to an internal 'C'
%   routine.
% 
% Examples:
% 
%   Increase the NURBS surface twice along the V direction.
%   esrf = nrbdegelev(srf, [0, 2]); 
% 
% See also:
% 
%   bspdegelev
%
%    Copyright (C) 2000 Mark Spink, 2010 Rafel Vazquez
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

if nargin < 2
  error('Input argument must include the NURBS and degree increment.');
end

if ~isstruct(nurbs)
  error('NURBS representation is not structure!');
end

if ~strcmp(nurbs.form,'B-NURBS')
  error('Not a recognised NURBS representation');
end

degree = nurbs.order-1;

if iscell(nurbs.knots)
 if size(nurbs.knots,2) == 3
  % NURBS represents a volume
  [dim,num1,num2,num3] = size(nurbs.coefs);

  % Degree elevate along the w direction
  if ntimes(3) == 0
    coefs = nurbs.coefs;
    knots{3} = nurbs.knots{3};
  else
    coefs = reshape(nurbs.coefs,4*num1*num2,num3);
    [coefs,knots{3}] = bspdegelev(degree(3),coefs,nurbs.knots{3},ntimes(3));
    num3 = size(coefs,2);
    coefs = reshape(coefs,[4 num1 num2 num3]);
  end

  % Degree elevate along the v direction
  if ntimes(2) == 0
    knots{2} = nurbs.knots{2};
  else
    coefs = permute(coefs,[1 2 4 3]);
    coefs = reshape(coefs,4*num1*num3,num2);
    [coefs,knots{2}] = bspdegelev(degree(2),coefs,nurbs.knots{2},ntimes(2));
    num2 = size(coefs,2);
    coefs = reshape(coefs,[4 num1 num3 num2]);
    coefs = permute(coefs,[1 2 4 3]);
  end

  % Degree elevate along the u direction
  if ntimes(1) == 0
    knots{1} = nurbs.knots{1};
  else
    coefs = permute(coefs,[1 3 4 2]);
    coefs = reshape(coefs,4*num2*num3,num1);
    [coefs,knots{1}] = bspdegelev(degree(1),coefs,nurbs.knots{1},ntimes(1));
    coefs = reshape(coefs,[4 num2 num3 size(coefs,2)]);
    coefs = permute(coefs,[1 4 2 3]);
  end 

 elseif size(nurbs.knots,2) == 2
  % NURBS represents a surface
  [dim,num1,num2] = size(nurbs.coefs);

  % Degree elevate along the v direction
  if ntimes(2) == 0
    coefs = nurbs.coefs;
    knots{2} = nurbs.knots{2};
  else
    coefs = reshape(nurbs.coefs,4*num1,num2);
    [coefs,knots{2}] = bspdegelev(degree(2),coefs,nurbs.knots{2},ntimes(2));
    num2 = size(coefs,2);
    coefs = reshape(coefs,[4 num1 num2]);
  end

  % Degree elevate along the u direction
  if ntimes(1) == 0
    knots{1} = nurbs.knots{1};
  else
    coefs = permute(coefs,[1 3 2]);
    coefs = reshape(coefs,4*num2,num1);
    [coefs,knots{1}] = bspdegelev(degree(1),coefs,nurbs.knots{1},ntimes(1));
    coefs = reshape(coefs,[4 num2 size(coefs,2)]);
    coefs = permute(coefs,[1 3 2]);
  end 
 end
else

  % NURBS represents a curve
  if isempty(ntimes)
    coefs = nurbs.coefs;
    knots = nurbs.knots;
  else
    [coefs,knots] = bspdegelev(degree,nurbs.coefs,nurbs.knots,ntimes);
  end
  
end

% construct new NURBS
inurbs = nrbmak(coefs,knots);

end

%!demo
%! crv = nrbtestcrv;
%! plot(crv.coefs(1,:),crv.coefs(2,:),'bo')
%! title('Degree elevation along test curve: curve and control polygons.');
%! hold on;
%! plot(crv.coefs(1,:),crv.coefs(2,:),'b--');
%! nrbplot(crv,48);
%!
%! icrv = nrbdegelev(crv, 1);
%!
%! plot(icrv.coefs(1,:),icrv.coefs(2,:),'ro')
%! plot(icrv.coefs(1,:),icrv.coefs(2,:),'r--');
%! 
%! hold off;
