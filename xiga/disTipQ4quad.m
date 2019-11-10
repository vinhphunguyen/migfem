function [W,Q]=disTipQ4quad(order,phi,nodes,tip)
% Integration points for Q4 element contains a crack tip
% Input:
%  - order: order of Gauss rule for each subtriangle
%  - phi  : level sets define the crack
%  - nodes: global coordinates of element nodes
%  - tip  : global coordinates of the crack tip
% 
% Vinh Phu Nguyen
% nvinhphu@gmail.com

epsilon = 0.00001;
corner  = [1 2 3 4 1];
node    = [-1 -1; 1 -1; 1 1; -1 1];

coord = zeros(1,2);
ksi   = 0;
eta   = 0;
iter  = 10;

xnode = nodes(:,1);
ynode = nodes(:,2);

inc = 1;
while (inc < iter)
    [N,dNdxi]=lagrange_basis('Q4',coord);   % compute shape functions

    x     = N'*xnode;
    y     = N'*ynode;
    df1dr = dNdxi(:,1)' * xnode;
    df1ds = dNdxi(:,2)' * xnode;
    df2dr = dNdxi(:,1)' * ynode;
    df2ds = dNdxi(:,2)' * ynode;
 
    f1 = x - tip(1);
    f2 = y - tip(2);

    detF   = df1dr*df2ds - df1ds*df2dr ;
    dFinv  = 1.0/detF;
    
    invf(1,1) =  dFinv * df2ds;
    invf(1,2) = -dFinv * df1ds;
    invf(2,1) = -dFinv * df2dr;
    invf(2,2) =  dFinv * df1dr;

    ksi       = ksi - invf(1,1)*f1 - invf(1,2)*f2;
    eta       = eta - invf(2,1)*f1 - invf(2,2)*f2;

    coord(1) = ksi;
    coord(2) = eta;

    if( (abs(f1) < epsilon) && (abs(f2) < epsilon) )
        inc  = iter + 1;
        ntip = coord;
    else
        inc = inc + 1;
    end
end

% loop on element edges
for i = 1:4
    n1 = corner(i);
    n2 = corner(i+1);
    if phi(n1)*phi(n2) < 0
        r    = phi(n1)/(phi(n1)-phi(n2));
        pnt  = (1-r)*node(n1,:)+r*node(n2,:);
        node = [node;pnt];
    end
end

% insert the tip into the Delaunay triangulation
node = [node;ntip];

% get decompused triangles
tri = delaunay(node(:,1),node(:,2));
tri = tricheck(node,tri);

% loop over subtriangles to get quadrature points and weights
pt = 1;
for e=1:size(tri,1)
    [w,q] = quadrature(order,'TRIANGULAR',2);
    % transform quadrature points into the parent element
    coord = node(tri(e,:),:);
    a = det([coord,[1;1;1]])/2;
    if ( a<0 )  % need to swap connectivity
        coord=[coord(2,:);coord(1,:);coord(3,:)];
        a = det([coord,[1;1;1]])/2;
    end

    if ( a~=0 )
        for n=1:length(w)
            N = lagrange_basis('T3',q(n,:));
            Q(pt,:) = N'*coord;
            W(pt,1) = 2*w(n)*a;
            pt = pt+1;
        end
    end
end








