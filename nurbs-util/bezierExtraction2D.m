function [C,Cxi,Cet] = bezierExtraction2D(uknot,vknot,p,q)
%
% Bezier extraction operators for a 2D NURBS/BSpline.
%
% VP Nguyen
% Cardiff University, UK

% Bezier extraction operators for xi and eta
% nb1: number of elements along xi direction
% nb2: number of elements along eta direction

[Cxi,nb1]  = bezierExtraction(uknot,p);
[Cet,nb2]  = bezierExtraction(vknot,q);

% For Bsplines/NURBS, the element Bezier extraction
% operator is square.

size1 = size(Cxi(:,:,1),1);
size2 = size(Cet(:,:,1),1);

% Bezier extraction operators for the whole mesh
% as the tensor product of Cxi and Cet

C = zeros(size1*size2,size1*size2,nb1*nb2);

for eta=1:nb2
    for xi=1:nb1
        e = (eta-1)*nb1 + xi;
        for row=1:size2
            ird = (row-1)*size1 + 1;
            jrd =  row*size1;
            for col=1:size2
                icd = (col-1)*size1 + 1;
                jcd =  col*size1;
                C(ird:jrd,icd:jcd,e) = Cet(row,col,eta) * Cxi(:,:,xi);
            end
        end
    end
end