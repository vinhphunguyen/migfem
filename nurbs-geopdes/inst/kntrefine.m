% KNTREFINE: Refine a given knot vector by dividing each interval uniformly,
%             maintaining the continuity in previously existing knots.
%
%   [rknots]                  = kntrefine (knots, n_sub, degree, regularity)
%   [rknots, zeta]            = kntrefine (knots, n_sub, degree, regularity)
%   [rknots, zeta, new_knots] = kntrefine (knots, n_sub, degree, regularity)
%
% INPUT:
%
%     knots:      initial knot vector.
%     n_sub:      number of new knots to be added in each interval.
%     degree:     polynomial degree of the refined knot vector
%     regularity: maximum global regularity 
%
% OUTPUT:
%
%     rknots:    refined knot vector
%     zeta:      refined knot vector without repetitions
%     new_knots: new knots, to apply the knot insertion
%
% The regularity at the new inserted knots is the one given by the user.
% At previously existing knots, the regularity is the minimum
%  between the previous regularity, and the one given by the user.
%  This ensures optimal convergence rates in the context of IGA.
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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

function varargout = kntrefine (knots, n_sub, degree, regularity)

  if (iscell(knots))
    if (numel(n_sub)~=numel(degree) || numel(n_sub)~=numel(regularity) || ...
        numel(n_sub)~=numel(knots))
      error('kntrefine: n_sub, degree and regularity must have the same length as the number of knot vectors')
    end
    aux_knots = knots;
  else
    if (numel(n_sub)~=numel(degree) || numel(n_sub)~=numel(regularity) || ...
        numel(n_sub)~=1)
      error('kntrefine: n_sub, degree and regularity must have the same length as the number of knot vectors')
    end
    aux_knots = {knots};
  end

  if (nargout == 3)
    for idim = 1:numel(n_sub)
      if (degree(idim)+1 ~= sum (aux_knots{idim}==aux_knots{idim}(1)))
        error ('kntrefine: new_knots is only computed when the degree is maintained');
      end
    end
    for idim = 1:numel(n_sub)
      min_mult     = degree(idim) - regularity(idim);
      z            = unique (aux_knots{idim});
      nz           = numel (z);
      deg          = sum (aux_knots{idim} == z(1)) - 1;
      rknots{idim} = z(ones(1, degree(idim)+1));
      new_knots{idim} = [];
 
      for ik = 2:nz
        insk = linspace (z(ik-1), z(ik), n_sub(idim) + 2);
        insk = vec (repmat (insk(2:end-1), min_mult, 1))';
        old_mult = sum (aux_knots{idim} == z(ik));
        mult = max (min_mult, degree(idim) - deg + old_mult);
        rknots{idim} = [rknots{idim}, insk, z(ik*ones(1, mult))];
        new_knots{idim} = [new_knots{idim}, insk, z(ik*ones(1, mult-old_mult))];
      end
      zeta{idim} = unique (rknots{idim});
    end
    if (~iscell(knots))
      rknots = rknots{1};
      zeta = zeta{1};
      new_knots = new_knots{1};
    end
    varargout{1} = rknots;
    varargout{2} = zeta;
    varargout{3} = new_knots;
  else
    for idim = 1:numel(n_sub)
      min_mult     = degree(idim) - regularity(idim);
      z            = unique (aux_knots{idim});
      nz           = numel (z);
      deg          = sum (aux_knots{idim} == z(1)) - 1;
      rknots{idim} = z(ones(1, degree(idim)+1));
 
      for ik = 2:nz
        insk = linspace (z(ik-1), z(ik), n_sub(idim) + 2);
        insk = vec (repmat (insk(2:end-1), min_mult, 1))';
        old_mult = sum (aux_knots{idim} == z(ik));
        mult = max (min_mult, degree(idim) - deg + old_mult);
        rknots{idim} = [rknots{idim}, insk, z(ik*ones(1, mult))];
      end
      zeta{idim} = unique (rknots{idim});
    end
    if (~iscell(knots))
      rknots = rknots{1};
      zeta = zeta{1};
    end
    varargout{1} = rknots;
    if (nargout == 2)
      varargout{2} = zeta;
    end
  end
end

function v = vec (in)
  v = in(:);
end

%!shared nrbs
%!test
%! knots = {[0 0 1 1] [0 0 0 1 1 1]};
%! coefs(1,:,:) = [1 sqrt(2)/2 0; 2 sqrt(2) 0];
%! coefs(2,:,:) = [0 sqrt(2)/2 1; 0 sqrt(2) 2];
%! coefs(4,:,:) = [1 sqrt(2)/2 1; 1 sqrt(2)/2 1];
%! nrbs = nrbmak (coefs, knots);
%! nrbs = nrbkntins (nrbs, {[] [0.5 0.6 0.6]});
%! nrbs = nrbdegelev (nrbs, [0 1]);
%! nrbs = nrbkntins (nrbs, {[] [0.4]});
%! rknots = kntrefine (nrbs.knots, [1 1], [1 1], [0 0]);
%! assert (rknots{1} == [0 0 0.5 1 1]);
%! assert (rknots{2} == [0 0 0.2 0.4 0.45 0.5 0.55 0.6 0.8 1 1]);
%!
%!test
%! rknots = kntrefine (nrbs.knots, [1 1], [3 3], [0 0]);
%! assert (rknots{1}, [0 0 0 0 0.5 0.5 0.5 1 1 1 1]);
%! assert (rknots{2}, [0 0 0 0 0.2 0.2 0.2 0.4 0.4 0.4 0.45 0.45 0.45 0.5 0.5 0.5 0.55 0.55 0.55 0.6 0.6 0.6 0.8 0.8 0.8 1 1 1 1]);
%!
%!test
%! rknots = kntrefine (nrbs.knots, [1 1], [3 3], [2 2]);
%! assert (rknots{1}, [0 0 0 0 0.5 1 1 1 1]);
%! assert (rknots{2}, [0 0 0 0 0.2 0.4 0.45 0.5 0.5 0.55 0.6 0.6 0.6 0.8 1 1 1 1]);
%!
%!test
%! rknots = kntrefine (nrbs.knots, [1 1], [4 4], [0 0]);
%! assert (rknots{1}, [0 0 0 0 0 0.5 0.5 0.5 0.5 1 1 1 1 1]);
%! assert (rknots{2}, [0 0 0 0 0 0.2 0.2 0.2 0.2 0.4 0.4 0.4 0.4 0.45 0.45 0.45 0.45 0.5 0.5 0.5 0.5 0.55 0.55 0.55 0.55 0.6 0.6 0.6 0.6 0.8 0.8 0.8 0.8 1 1 1 1 1]);
%!
%!test
%! rknots = kntrefine (nrbs.knots, [1 1], [4 4], [3 3]);
%! assert (rknots{1}, [0 0 0 0 0 0.5 1 1 1 1 1]);
%! assert (rknots{2}, [0 0 0 0 0 0.2 0.4 0.4 0.45 0.5 0.5 0.5 0.55 0.6 0.6 0.6 0.6 0.8 1 1 1 1 1]);
%!
%!test
%! knots = [0 0 0 0 0.4 0.5 0.5 0.6 0.6 0.6 1 1 1 1];
%! rknots = kntrefine (knots, 1, 4, 3);
%! assert (rknots, [0 0 0 0 0 0.2 0.4 0.4 0.45 0.5 0.5 0.5 0.55 0.6 0.6 0.6 0.6 0.8 1 1 1 1 1]);
