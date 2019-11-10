function srf = nrbextrude(curve,vector)

%
% NRBEXTRUDE: Construct a NURBS surface by extruding a NURBS curve, or 
%  construct a NURBS volume by extruding a NURBS surface.
% 
% Calling Sequence:
% 
%   srf = nrbextrude(crv,vec);
% 
% INPUT:
% 
%   crv		: NURBS curve or surface to extrude, see nrbmak.
% 
%   vec		: Vector along which the entity is extruded.
%
% OUTPUT: 
% 
%   srf		: NURBS surface or volume constructed.
% 
% Description:
% 
%   Constructs either a NURBS surface by extruding a NURBS curve along a  
%   defined vector, or a NURBS volume by extruding a NURBS surface. In the 
%   first case, the NURBS curve forms the U direction of the surface edge, and
%   is extruded along the vector in the V direction. In the second case, the 
%   original surface forms the U and V direction of the volume, and is extruded
%   along the W direction.
%
% Examples:
% 
%   Form a hollow cylinder by extruding a circle along the z-axis.
%
%   srf = nrbextrude(nrbcirc, [0,0,1]);
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

if (nargin < 2)
  error('Error too few input arguments!');
end

if (iscell (curve.knots))
  if (numel (curve.knots) == 3)
    error('Nurbs volumes cannot be extruded!');
  end
  for ii = 1:size(curve.coefs,3)
    coefs(:,:,ii) = vectrans(vector) * squeeze (curve.coefs(:,:,ii));
  end
  coefs = cat(4,curve.coefs,coefs);
  srf = nrbmak(coefs,{curve.knots{:}, [0 0 1 1]});
else
  coefs = cat(3,curve.coefs,vectrans(vector)*curve.coefs);
  srf = nrbmak(coefs,{curve.knots, [0 0 1 1]});
end

end

%!demo
%! crv = nrbtestcrv;
%! srf = nrbextrude(crv,[0 0 5]);
%! nrbplot(srf,[40 10]);
%! title('Extrusion of a test curve along the z-axis');
%! hold off
%
%!demo
%! crv1 = nrbcirc (1, [0 0], 0, pi/2);
%! crv2 = nrbcirc (2, [0 0], 0, pi/2);
%! srf  = nrbruled (crv1, crv2);
%! vol  = nrbextrude (srf, [0 0 1]);
%! nrbplot (vol, [30 10 10])
%! title ('Extrusion of the quarter of a ring')
%
%!demo
%! srf = nrbtestsrf;
%! vol = nrbextrude(srf, [0 0 10]);
%! nrbplot(vol,[20 20 20]);
%! title('Extrusion of a test surface along the z-axis');
%! hold off
