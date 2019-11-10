function inurbs = nrbkntins(nurbs,iknots)
% 
% NRBKNTINS: Insert a single or multiple knots into a NURBS curve,
%            surface or volume.
% 
% Calling Sequence:
% 
%   icrv = nrbkntins(crv,iuknots);
%   isrf = nrbkntins(srf,{iuknots ivknots});
%   ivol = nrbkntins(vol,{iuknots ivknots iwknots});
% 
% INPUT:
% 
%   crv		: NURBS curve, see nrbmak.
% 
%   srf		: NURBS surface, see nrbmak.
% 
%   srf		: NURBS volume, see nrbmak.
% 
%   iuknots	: Knots to be inserted along U direction.
% 
%   ivknots	: Knots to be inserted along V direction.
% 
%   iwknots	: Knots to be inserted along W direction.
% 
% OUTPUT:
% 
%   icrv	: new NURBS structure for a curve with knots inserted.
% 
%   isrf	: new NURBS structure for a surface with knots inserted.
% 
%   ivol	: new NURBS structure for a volume with knots inserted.
% 
% Description:
% 
%   Inserts knots into the NURBS data structure, these can be knots at
%   new positions or at the location of existing knots to increase the
%   multiplicity. Note that the knot multiplicity cannot be increased
%   beyond the order of the spline. Knots along the V direction can only
%   inserted into NURBS surfaces, not curve that are always defined along
%   the U direction. This function use the B-Spline function bspkntins,
%   which interfaces to an internal 'C' routine.
% 
% Examples:
% 
%   Insert two knots into a curve, one at 0.3 and another
%   twice at 0.4
%
%   icrv = nrbkntins(crv, [0.3 0.4 0.4])
% 
%   Insert into a surface two knots as (1) into the U knot
%   sequence and one knot into the V knot sequence at 0.5.
%
%   isrf = nrbkntins(srf, {[0.3 0.4 0.4] [0.5]})
% 
% See also:
% 
%   bspkntins
%
% Note:
%
%   No knot multiplicity will be increased beyond the order of the spline.
%
%    Copyright (C) 2000 Mark Spink, 2010 Rafael Vazquez
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
  error('Input argument must include the NURBS and knots to be inserted');
end

if ~isstruct(nurbs)
  error('NURBS representation is not structure!');
end

if ~strcmp(nurbs.form,'B-NURBS')
  error('Not a recognised NURBS representation');
end

degree = nurbs.order-1;

if iscell(nurbs.knots)
 if size(nurbs.knots,2)==3
  % NURBS represents a volume
  num1 = nurbs.number(1);
  num2 = nurbs.number(2);
  num3 = nurbs.number(3);

  % Insert knots along the w direction
  if isempty(iknots{3})
    coefs = nurbs.coefs;
    knots{3} = nurbs.knots{3};
  else
    coefs = reshape(nurbs.coefs,4*num1*num2,num3);
    [coefs,knots{3}] = bspkntins(degree(3),coefs,nurbs.knots{3},iknots{3});
    num3 = size(coefs,2);
    coefs = reshape(coefs,[4 num1 num2 num3]);
  end

  % Insert knots along the v direction
  if isempty(iknots{2})
    knots{2} = nurbs.knots{2};
  else
    coefs = permute(coefs,[1 2 4 3]);
    coefs = reshape(coefs,4*num1*num3,num2);
    [coefs,knots{2}] = bspkntins(degree(2),coefs,nurbs.knots{2},iknots{2});
    num2 = size(coefs,2);
    coefs = reshape(coefs,[4 num1 num3 num2]);
    coefs = permute(coefs,[1 2 4 3]);
  end

  % Insert knots along the u direction
  if isempty(iknots{1})
    knots{1} = nurbs.knots{1};
  else   
    coefs = permute(coefs,[1 3 4 2]);
    coefs = reshape(coefs,4*num2*num3,num1);
    [coefs,knots{1}] = bspkntins(degree(1),coefs,nurbs.knots{1},iknots{1});
    coefs = reshape(coefs,[4 num2 num3 size(coefs,2)]);
    coefs = permute(coefs,[1 4 2 3]);
  end
     
 elseif size(nurbs.knots,2)==2
  % NURBS represents a surface
  num1 = nurbs.number(1);
  num2 = nurbs.number(2);

  % Insert knots along the v direction
  if isempty(iknots{2})
    coefs = nurbs.coefs;
    knots{2} = nurbs.knots{2};
  else
    coefs = reshape(nurbs.coefs,4*num1,num2);
    [coefs,knots{2}] = bspkntins(degree(2),coefs,nurbs.knots{2},iknots{2});
    num2 = size(coefs,2);
    coefs = reshape(coefs,[4 num1 num2]);
  end

  % Insert knots along the u direction
  if isempty(iknots{1})
    knots{1} = nurbs.knots{1};
  else   
    coefs = permute(coefs,[1 3 2]);
    coefs = reshape(coefs,4*num2,num1);
    [coefs,knots{1}] = bspkntins(degree(1),coefs,nurbs.knots{1},iknots{1});
    coefs = reshape(coefs,[4 num2 size(coefs,2)]);
    coefs = permute(coefs,[1 3 2]);
  end
 end
else

  % NURBS represents a curve
  if isempty(iknots)
    coefs = nurbs.coefs;
    knots = nurbs.knots;
  else
    [coefs,knots] = bspkntins(degree,nurbs.coefs,nurbs.knots,iknots);  
  end

end

% construct new NURBS
inurbs = nrbmak(coefs,knots); 

end

%!demo
%! crv = nrbtestcrv;
%! plot(crv.coefs(1,:),crv.coefs(2,:),'bo')
%! title('Knot insertion along test curve: curve and control polygons.');
%! hold on;
%! plot(crv.coefs(1,:),crv.coefs(2,:),'b--');
%!
%! nrbplot(crv,48);
%!
%! icrv = nrbkntins(crv,[0.125 0.375 0.625 0.875] );
%! plot(icrv.coefs(1,:),icrv.coefs(2,:),'ro')
%! plot(icrv.coefs(1,:),icrv.coefs(2,:),'r--');
%! hold off