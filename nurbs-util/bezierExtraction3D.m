function [C,Cxi,Cet,Cze] = bezierExtraction3D(uKnot,vKnot,wKnot,p,q,r)
%
% Bezier extraction operators for a 3D NURBS/BSpline.
% Using [C nb] = bezierExtraction(knot,p) for each direction
% and tensor product Ci * Cj * Ck for the final trivariate C.
%
% VP Nguyen
% Cardiff University, UK

% Bezier extraction operators for xi x eta

[Cxe,Cxi,Cet]  = bezierExtraction2D(uKnot,vKnot,p,q);

% Bezier extractor for the third direction

[Cze,nb3]      = bezierExtraction(wKnot,r);

% For Bsplines/NURBS, the element Bezier extraction
% operator is square.

size12 = size(Cxi(:,:,1),1) * size(Cet(:,:,1),1);
size3  = size(Cze(:,:,1),1);
nb12   = size(Cxe,3);

% Bezier extraction operators for the whole mesh
% as the tensor product of Cxi and Cet

C = zeros(size12*size3,size12*size3,nb12*nb3);


for zeta=1:nb3
    for xi=1:nb12
        e = (zeta-1)*nb12 + xi;
        for row=1:size3
            ird = (row-1)*size12 + 1;
            jrd =  row*size12;
            for col=1:size3
                icd = (col-1)*size12 + 1;
                jcd =  col*size12;
                C(ird:jrd,icd:jcd,e) = Cze(row,col,zeta) * Cxe(:,:,xi);
            end
        end
    end
end