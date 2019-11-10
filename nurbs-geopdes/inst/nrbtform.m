function nurbs = nrbtform(nurbs,tmat)
% 
% NRBTFORM: Apply transformation matrix to the NURBS.
% 
% Calling Sequence:
% 
%   tnurbs = nrbtform(nurbs,tmatrix);
% 
% INPUT:
% 
%   nurbs	: NURBS data structure (see nrbmak for details).
% 
%   tmatrix     : Transformation matrix, a matrix of size (4,4) defining
%                 a single or multiple transformations.
%
% OUTPUT:
%
%   tnurbs	: The return transformed NURBS data structure.
% 
% Description:
% 
%   The NURBS is transform as defined a transformation matrix of size (4,4),
%   such as a rotation, translation or change in scale. The transformation
%   matrix can define a single transformation or multiple series of
%   transformations. The matrix can be simple constructed by the functions
%   vecscale, vectrans, vecrotx, vecroty, and vecrotz.
%     
% Examples:
% 
%   Rotate a square by 45 degrees about the z axis.
%
%   rsqr = nrbtform(nrbrect(), vecrotz(deg2rad(45)));
%   nrbplot(rsqr, 1000);
% 
% See also:
% 
%   vecscale, vectrans, vecrotx, vecroty, vecrotz
%
%    Copyright (C) 2000 Mark Spink
%    Copyright (C) 2010 Rafael Vazquez
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
  error('Not enough input arguments!');
end;

if iscell(nurbs.knots)
 if size(nurbs.knots,2) == 2
  % NURBS is a surface
  [dim,nu,nv] = size(nurbs.coefs);
  nurbs.coefs = reshape(tmat*reshape(nurbs.coefs,dim,nu*nv),[dim nu nv]);
 elseif size(nurbs.knots,2) == 3
  % NURBS is a volume
  [dim,nu,nv,nw] = size(nurbs.coefs);
  nurbs.coefs = reshape(tmat*reshape(nurbs.coefs,dim,nu*nv*nw),[dim nu nv nw]);
 end
else
  % NURBS is a curve
  nurbs.coefs = tmat*nurbs.coefs;
end

end

%!demo
%! xx = vectrans([2.0 1.0])*vecroty(pi/8)*vecrotx(pi/4)*vecscale([1.0 2.0]);
%! c0 = nrbtform(nrbcirc, xx);
%! nrbplot(c0,50);
%! grid on
%! title('Construction of an ellipse by transforming a unit circle.');
%! hold off