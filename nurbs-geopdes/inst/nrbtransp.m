function tsrf = nrbtransp(srf)
% 
% NRBTRANSP: Transpose a NURBS surface, by swapping U and V directions.
% 
% Calling Sequence:
% 
%   tsrf = nrbtransp(srf)
%
% INPUT:
% 
%   srf		: NURBS surface, see nrbmak.
%
% OUTPUT:
% 
%   tsrf	: NURBS surface with U and V diretions transposed.
% 
% Description:
% 
%   Utility function that transposes a NURBS surface, by swapping U and
%   V directions. NURBS curves cannot be transposed.
%
%    Copyright (C) 2000 Mark Spink
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

if ~iscell(srf.knots)
  error(' A NURBS curve cannot be transposed.');
elseif size(srf.knots,2) == 3
  error('The transposition of NURBS volumes has not been implemented.');
end  

tsrf = nrbmak(permute(srf.coefs,[1 3 2]), fliplr(srf.knots));

end

%!demo
%! srf = nrb4surf([0 0 0], [1 0 1], [0 1 1], [1 1 2]);
%! nrbplot(srf,[20 5]);
%! title('Plane surface and its transposed (translated)')
%! hold on
%! srf.coefs(3,:,:) = srf.coefs(3,:,:) + 10;
%! srf = nrbtransp(srf);
%! nrbplot(srf,[20 5]);
%! hold off