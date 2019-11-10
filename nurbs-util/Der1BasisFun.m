function ders = Der1BasisFun(i,u,p,U)
%--------------------------------------------------------------
%function ders = Der1BasisFun(i,u,p,U)
% NURBS-Book modified (algorithm A2.3)
% evalute nonzero basis functions and first derivative
%INPUT:
% i          : current knotspan
% u          : evaluation point
% p          : degree of the basis functions
% U          : knot vector
%OUTPUT:
% ders       : matrix (2 , p+1)
%              1. row vector (dim p+1)
%              values of the basis function N_(i-p) ... N_(i)
%              at the evaluation point u
%              2. first derivation
%--------------------------------------------------------------

ders   = zeros(2,p+1);
N      = zeros(p+1,p+1);
N(1,1) = 1;
left   = zeros(1,p+1);
right  = zeros(1,p+1);

for j=1:p
    left(j+1)  = u-U(i+1-j+1);
    right(j+1) = U(i+j+1)-u;
    saved = 0;
    for r=0:j-1
        N(j+1,r+1) = right(r+2) + left(j-r+1);
        temp = N(r+1,j)/N(j+1,r+1);
        N(r+1,j+1) = saved + right(r+2)*temp;
        saved = left(j-r+1)*temp;
    end
    N(j+1,j+1) = saved;
end

for j=0:p
    ders(1,j+1) = N(j+1,p+1);
end
%derivative
for r=0:p
    %kth derivative
    if(r>=1)
        ders(2,r+1)= N(r,p)/N(p+1,r);
    end
    if (r<= p-1)
        ders(2,r+1) = ders(2,r+1) -N(r+1,p)/N(p+1,r+1);
    end
end
% Multiply through by the correct factors
ders(2,:) = ders(2,:)*p;
end