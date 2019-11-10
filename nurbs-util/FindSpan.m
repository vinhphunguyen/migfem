function knotSpanIndex = FindSpan(n,p,u,U)
%--------------------------------------------------------------
%function knotSpanIndex = FindSpan(n,p,u,U)
% NURBS-Book (algorithm A2.1)
% find the knot span index for one variable u
%INPUT:
% n          : number of basis function -1
%        NURBS-Book: np+1 # basis, np max index (startindex 0)
%        here        np   # basis and max index (startindex 1)
% p          : degree of the basis functions
% u          : evaluation point
% U          : knot vector (row vector)
%OUTPUT:
% knotSpanIndex : index of knot span
%--------------------------------------------------------------


if (u == U(n+2))
    knotSpanIndex= n;
    return
end
low = p;
high = n+1;
mid = floor((low + high)/2);
while (u <U(mid+1) || u >= U(mid+2) )
    if( u < U(mid+1))
        high = mid;
    else
        low = mid;
    end
    mid = floor((low+high)/2);
end
knotSpanIndex = mid;
end