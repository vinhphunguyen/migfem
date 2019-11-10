function B = getBmatrix2D (dRdx)

nn2 = 2*size(dRdx,2);

B(1,1:2:nn2)  = dRdx(1,:);
B(2,2:2:nn2)  = dRdx(2,:);
B(3,1:2:nn2)  = dRdx(2,:);
B(3,2:2:nn2)  = dRdx(1,:);