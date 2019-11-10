function [c,h]=contoura(x,y,z,a,b),
% CONTOURA CONTOURA plot.
%       CONTOURA(Z) is a contour plot of matrix Z treating the values in Z
%       as heights above a plane.
%       CONTOURA(X,Y,Z), where X and Y are vectors and Z is a matrix, specifies
%       the X- and Y-axes used on the plot.
%       CONTOURA(X,Y,Z), where X, Y, and Z are vectors, interpolates
%       the data onto a rectangular grid.
%       CONTOURA(Z,N) and CONTOURA(X,Y,Z,N) draw N contour lines, 
%       overriding the default automatic value.
%       CONTOURA(Z,V) and CONTOURA(X,Y,Z,V) draw LENGTH(V) contour lines 
%       at the values specified in vector V.
%       CONTOURA(...,'linetype') draws with the color and linetype specified,
%       as in the PLOT command.
%
%       C = CONTOURA(...) returns contour matrix C as described in CONTOURC
%       and used by CLABEL.
%       [C,H] = CONTOURA(...) returns a column vector H of handles to LINE
%       objects, one handle per line. 
% 
%       NOTE: Do not rename this CONTOUR, it calls CONTOUR.
%       See also CLABEL, CONTOURC, CONTOUR3, GRADIENT, QUIVER, PRISM.

%       D. Thomas 11/17/94
%       Copyright (c) 1984-94 by The MathWorks, Inc.
A = 'xyzab';
C = ',,,, ';

args = [A(1:nargin);[C(1:nargin-1) ' ']];
args = [args(:)' ' '];

if (nargout == 1),
        argo = 'c=';
elseif (nargout == 2),
        argo = '[c,h]=';
else
        argo = '';
end 

if (nargin > 2)
        if (min(size(x)) == 1 & min(size(y)) == 1 & min(size(z)) == 1),
                
                if  ~(size(x(:),1) == size(y(:),1) & size(y(:),1) == size(z(:),1)),
                        error('Sizes of X, Y, and Z must be equal');
                end
                xmin=min(x(:));
                ymin=min(y(:));
                xmax=max(x(:));
                ymax=max(y(:));
                N = ceil(length(x(:))^(1/2));
                Xl = xmin:((xmax-xmin)/(N-1)):xmax;
                Yl = ymin:((ymax-ymin)/(N-1)):ymax;
                [Xi,Yi]=meshgrid(Xl,Yl);
                Z = griddata(x,y,z,Xi,Yi);
                if (length(args) > 7),
                        args = [',',args];
                end
                eval([argo,'contour(Xl,Yl,Z',args(7:length(args)),');']);
        else,
                eval([argo,'contour(' args ');']);
        end
else,
        eval([argo,'contour(' args ');']);
end



