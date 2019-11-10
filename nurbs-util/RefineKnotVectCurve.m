function [Ubar,Qw] = RefineKnotVectCurve(n,p,U,Pw,X,r)
%--------------------------------------------------------------
%function [Ubar,Qw] = RefineKnotVectCurve(n,p,U,Pw,X,r)
% NURBS-Book (algorithm A5.4) (modified)
% insert multiple knots into curve
%INPUT:
% n         : number ob basis functions -1 !
%        NURBS-Book: n+1 # basis, np max index (startindex 0)
%        here        n   # basis and max index (startindex 1)
% p          : degree of the basis functions
% U         : old knotvector
% Pw         : old control points
% X          : vector of new knots (multiple entries possible)
% r          :  (size of X) -1 (count the multple entries as well
%             reason: same as above: max X index
%OUTPUT:
% Ubar        : newknot vector
% Qw         : new control points
%--------------------------------------------------------------

%initialise arrays;
dim  = size(Pw,2);
Qw   = zeros(n+r+2,dim);
Ubar = zeros(1,n+p+1+r);
%
m = n+p+1;
a = FindSpan(n,p,X(1),U);
b = FindSpan(n,p,X(r+1),U);
b = b+1;

for j=0:a-p
    Qw(j+1,:) = Pw(j+1,:);
end
for j=b-1:n
    Qw(j+r+2,:) = Pw(j+1,:);
end
for j=0:a
    Ubar(j+1)= U(j+1);
end
for j=b+p :m
    Ubar(j+r+2) = U(j+1);
end
i=b+p-1;
k=b+p+r;
for j=r:-1:0
    while (X(j+1) <= U(i+1) && i>a)
        Qw(k-p,:) = Pw(i-p,:);
        Ubar(k+1) = U(i+1);
        k=k-1;
        i=i-1;
    end
    Qw(k-p,:) = Qw(k-p+1,:);
    for l=1:p
        ind = k-p+l;
        alfa = Ubar(k+l+1) - X(j+1);
        if (abs(alfa) == 0)
            Qw(ind,:) = Qw(ind+1,:);
        else
            alfa = alfa/(Ubar(k+l+1) - U(i-p+l+1));
            Qw(ind,:) = alfa* Qw(ind,:) + (1-alfa)* Qw(ind+1,:);
        end
    end
    Ubar(k+1) = X(j+1);
    k=k-1;
end
end