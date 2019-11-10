function S = SurfacePoint(n,p,U,m,q,V,P,dim,u,v)
%--------------------------------------------------------------
%function S = SurfacePoint(n,p,U,m,q,V,P,u,v)
% NURBS-Book (algorithm A4.5) (modified)
% calculate point on a surface
% this can be used for B-Spline surfaces
% or NURBS surfaces in projective coordinates (correct variable dim!)
%INPUT:
% n         : number ob basis functions -1 !  - x-direction
%        NURBS-Book: n+1 # basis, np max index (startindex 0)
%        here        n   # basis and max index (startindex 1)
% p          : degree of the basis functions - x-direction
% U          : knotvector - x-direction
% m          : number ob basis functions -1 !  - y-direction
% q          : degree of the basis functions - y-direction
% V          : knotvector - x-direction
% P          : control points
% dim        : dimension of control points
% u          : x-coordinate
% v          : y-coordinate
%OUTPUT:
% S          : coordinates of the point on the surface
%--------------------------------------------------------------

uspan = FindSpan(n,p,u,U);
vspan = FindSpan(m,q,v,V);

Nu    = BasisFun(uspan, u,p,U);
Nv    = BasisFun(vspan, v,q,V);

uind = uspan -p;
S    = zeros(1,dim);

for l=0:q
    temp = zeros(1,dim);
    vind = vspan-q+l;
    for i=0:p
        %access control Point P(uind+i+1,vind+1)
        CP   = P(uind+i+1 + vind*(n+1),:);
        temp = temp + Nu(i+1)*CP;
    end
    S = S + Nv(l+1)*temp;
end


