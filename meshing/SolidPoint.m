function S = SolidPoint(n,p,U,m,q,V,l,r,W,P,dim,u,v,w)
%--------------------------------------------------------------
%function S = SolidPoint(n,p,U,m,q,V,P,u,v)
% this can be used for B-Spline solid
% or NURBS solid in projective coordinates (correct variable dim!)
%INPUT:
% n         : number ob basis functions -1 !  - x-direction
%        NURBS-Book: n+1 # basis, np max index (startindex 0)
%        here        n   # basis and max index (startindex 1)
% p          : degree of the basis functions - x-direction
% U          : knotvector - x-direction
% m          : number ob basis functions -1 !  - y-direction
% q          : degree of the basis functions - y-direction
% V          : knotvector - y-direction
% P          : control points
% dim        : dimension of control points
% u          : xi-coordinate
% v          : eta-coordinate
% w          : zeta-coordinate
%OUTPUT:
% S          : coordinates of the point on the solid
%--------------------------------------------------------------

% find spans of (u,v,w)

uspan = FindSpan(n,p,u,U);
vspan = FindSpan(m,q,v,V);
wspan = FindSpan(l,r,w,W);

% compute non-zero B-spline basis functions

Nu    = BasisFun(uspan, u,p,U);
Nv    = BasisFun(vspan, v,q,V);
Nw    = BasisFun(wspan, w,r,W);

% compute point on solid using B-spline interpolation

uind = uspan -p;
S    = zeros(1,dim);

for k=0:r
  wind = wspan-r+k;
  for j=0:q
      vind = vspan-q+j;
      for i=0:p
	      CP   = P(uind+i+1 + (n+1) * ( (m+1)*wind + vind),:); 
          S    = S + Nu(i+1) * Nv(j+1) * Nw(k+1) * CP;
      end
  end
end


