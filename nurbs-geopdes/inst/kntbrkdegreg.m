% KNTBRKDEGREG: Construct an open knot vector by giving the sequence of
%                knots, the degree and the regularity.
%
%   knots = kntbrkdegreg (breaks, degree)
%   knots = kntbrkdegreg (breaks, degree, regularity)
%
% INPUT:
%
%     breaks:     sequence of knots.
%     degree:     polynomial degree of the splines associated to the knot vector.
%     regularity: splines regularity.
%
% OUTPUT:
%
%     knots:  knot vector.
%
% If REGULARITY has as many entries as BREAKS, or as the number of interior
%   knots, a different regularity will be assigned to each knot. If
%   REGULARITY is not present, it will be taken equal to DEGREE-1.
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

function knots = kntbrkdegreg (breaks, degree, reg)

if (iscell (breaks))
  if (nargin == 2)
    reg = degree - 1;
  end
  
  if (numel(breaks)~=numel(degree) || numel(breaks)~=numel(reg))
    error('kntbrkdegreg: degree and regularity must have the same length as the number of knot vectors')
  end

  degree = num2cell (degree);

  if (~iscell (reg))
    reg = num2cell (reg);
  end

  knots = cellfun (@do_kntbrkdegreg, breaks, degree, reg, 'uniformoutput', false);
else

  if (nargin == 2)
    reg = degree - 1;
  end

  knots = do_kntbrkdegreg (breaks, degree, reg);
end

end

function knots = do_kntbrkdegreg (breaks, degree, reg)

if (numel (breaks) < 2)
  error ('kntbrkdegreg: the knots sequence should contain at least two points') 
end

if (numel (reg) == 1)
  mults = [-1, (degree (ones (1, numel (breaks) - 2)) - reg), -1];
elseif (numel (reg) == numel (breaks))
  mults = degree - reg;
elseif (numel (reg) == numel (breaks) - 2)
  mults = [-1 degree-reg -1];
else
  error('kntbrkdegreg: the length of mult should be equal to one or the number of knots')
end

if (any (reg < -1))
  warning ('kntbrkdegreg: for some knots the regularity is lower than -1')
elseif (any (reg > degree-1))
  error('kntbrkdegreg: the regularity should be lower than the degree')
end

knots = kntbrkdegmult (breaks, degree, mults);

end

%!test
%! breaks = [0 1 2 3 4];
%! degree = 3;
%! knots = kntbrkdegreg (breaks, degree);
%! assert (knots, [0 0 0 0 1 2 3 4 4 4 4])

%!test
%! breaks = [0 1 2 3 4];
%! degree = 3;
%! reg    = 1;
%! knots = kntbrkdegreg (breaks, degree, reg);
%! assert (knots, [0 0 0 0 1 1 2 2 3 3 4 4 4 4])

%!test
%! breaks = [0 1 2 3 4];
%! degree = 3;
%! reg    = [0 1 2];
%! knots = kntbrkdegreg (breaks, degree, reg);
%! assert (knots, [0 0 0 0 1 1 1 2 2 3 4 4 4 4])

%!test
%! breaks = {[0 1 2 3 4] [0 1 2 3]};
%! degree = [3 2];
%! reg    = {[0 1 2] 0};
%! knots = kntbrkdegreg (breaks, degree, reg);
%! assert (knots, {[0 0 0 0 1 1 1 2 2 3 4 4 4 4] [0 0 0 1 1 2 2 3 3 3]})

