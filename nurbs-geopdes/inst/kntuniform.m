% KNTUNIFORM: generate uniform open knot vectors in the reference domain.
%
%   [csi, zeta] = kntuniform (num, degree, regularity)
%
% INPUT:
%     
%     num:        number of breaks (in each direction)
%     degree:     polynomial degree (in each direction)
%     regularity: global regularity (in each direction)
%
% OUTPUT:
%
%     csi:  knots
%     zeta: breaks = knots without repetitions
% 
% Copyright (C) 2009, 2010 Carlo de Falco
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

function [csi, zeta] = kntuniform (num, degree, regularity)
  
  if (numel(num)~=numel(degree) || numel(num)~=numel(regularity))
    error('kntuniform: num, degree and regularity must have the same length')
  else
    for idim=1:numel(num)
      zeta{idim} = linspace (0, 1, num(idim));
      rep  = degree(idim) - regularity(idim);
      if (rep > 0)
        csi{idim}  = [zeros(1, degree(idim)+1-rep)...
          reshape(repmat(zeta{idim}, rep, 1), 1, []) ones(1, degree(idim)+1-rep)];
      else
        error ('kntuniform: regularity requested is too high')
      end
    end
  end
end
