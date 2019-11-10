function rnrb = nrbreverse(nrb)
%
% NRBREVERSE: Reverse the evaluation direction of a NURBS curve or surface.
% 
% Calling Sequence:
% 
%   rnrb = nrbreverse(nrb);
% 
% INPUT:
% 
%   nrb		: NURBS data structure, see nrbmak.
%
% OUTPUT:
% 
%   rnrb	: Reversed NURBS.
% 
% Description:
% 
%   Utility function to reverse the evaluation direction of a NURBS
%   curve or surface.
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

if nargin ~= 1
  error('Incorrect number of input arguments');
end

if iscell(nrb.knots)
 if size(nrb.knots,2) == 3
  error('The function nrbreverse is not yet ready for volumes')
 else
  % reverse a NURBS surface
  coefs = nrb.coefs(:,:,end:-1:1);
  rnrb = nrbmak(coefs(:,end:-1:1,:), {1.0-fliplr(nrb.knots{1}),...
                1.0-fliplr(nrb.knots{2})});
 end

else

  % reverse a NURBS curve
  rnrb = nrbmak(fliplr(nrb.coefs), 1.0-fliplr(nrb.knots));

end

end

%!demo
%! pnts = [0.5 1.5 3.0 7.5 8.5;
%!         3.0 5.5 1.5 4.0 4.5;
%!         0.0 0.0 0.0 0.0 0.0];
%! crv1 = nrbmak(pnts,[0 0 0 1/2 3/4 1 1 1]);
%! crv2 = nrbreverse(crv1);
%! fprintf('Knots of the original curve\n')
%! disp(crv1.knots)
%! fprintf('Knots of the reversed curve\n')
%! disp(crv2.knots)
%! fprintf('Control points of the original curve\n')
%! disp(crv1.coefs(1:2,:))
%! fprintf('Control points of the reversed curve\n')
%! disp(crv2.coefs(1:2,:))
%! nrbplot(crv1,100)
%! hold on
%! nrbplot(crv2,100)
%! title('The curve and its reverse are the same')
%! hold off