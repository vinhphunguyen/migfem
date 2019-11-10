%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = BasisFun(i,u,p,U)
%--------------------------------------------------------------
%function N = BasisFun(i,p,u,U)
% NURBS-Book (algorithm A2.2)
% evalute nonzero basis functions
%INPUT:
% i          : current knotspan
% u          : evaluation point
% p          : degree of the basis functions
% U          : knot vector (row vector)
%OUTPUT:
% N          : row vector (dim p+1)
%              values of the basis function N_(i-p) ... N_(i)
%              at the evaluation point u
%--------------------------------------------------------------
N=zeros(1,p+1);
N(1)=1;
left=zeros(1,p+1);
right=zeros(1,p+1);
for j=1:p
    left(j+1) = u-U(i+1-j+1);
    right(j+1) = U(i+j+1)-u;
    saved = 0;
    for r=0:j-1
        temp = N(r+1)/(right(r+2)+left(j-r+1));
        N(r+1) = saved + right(r+2)*temp;
        saved = left(j-r+1)*temp;
    end
    N(j+1) = saved;
end
end