function nrbexport (nurbs, filename)

%
% NRBEXPORT: export NURBS geometries to a format compatible with the one used in GeoPDEs (version 0.6).
% 
% Calling Sequence:
% 
%   nrbexport (nurbs, filename);
% 
% INPUT:
% 
%   nurbs    : NURBS curve, surface or volume, see nrbmak.
%   filename : name of the output file.
% 
% 
% Description:
% 
%   The data of the nurbs structure is written in the file, in a format 
%     that can be read by GeoPDEs.
%
%    Copyright (C) 2011 Rafael Vazquez
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

fid = fopen (filename, 'w');
if (fid < 0)
  error ('nrbexport: cannot open file %s', filename);
end

ndim = numel (nurbs(1).order);
npatch = numel (nurbs);
fprintf (fid, '%s\n', '# nurbs mesh v.0.7');
fprintf (fid, '%s\n', '#');
fprintf (fid, '%s\n', ['# ' date]);
fprintf (fid, '%s\n', '#');

fprintf (fid, '%2i', ndim, 1);
fprintf (fid, '\n');
for iptc = 1:npatch
  fprintf (fid, '%s %i', 'PATCH', iptc);
  fprintf (fid, '\n');
  fprintf (fid, '%2i', nurbs(iptc).order-1);
  fprintf (fid, '\n');
  fprintf (fid, '%2i', nurbs(iptc).number);
  fprintf (fid, '\n');
  for ii = 1:ndim
    fprintf (fid, '%1.7f   ', nurbs(iptc).knots{ii});
    fprintf (fid, '\n');
  end

  if (ndim == 2)
    for ii = 1:ndim
      fprintf (fid, '%1.15f   ', nurbs(iptc).coefs(ii,:,:));
      fprintf (fid, '\n');
    end
    fprintf (fid, '%1.15f   ', nurbs(iptc).coefs(4,:,:));
    fprintf (fid, '\n');
  elseif (ndim == 3)
    for ii = 1:ndim
      fprintf (fid, '%1.15f   ', nurbs(iptc).coefs(ii,:,:,:));
      fprintf (fid, '\n');
    end
    fprintf (fid, '%1.15f   ', nurbs(iptc).coefs(4,:,:,:));
    fprintf (fid, '\n');
  end
end

end
