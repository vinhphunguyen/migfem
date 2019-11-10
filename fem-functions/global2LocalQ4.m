function [lpoint]=global2LocalQ4(nodes,gpoint)

epsilon = 0.00001;

coord = zeros(1,2);
ksi   = 0;
eta   = 0;
iter  = 10;

inc = 1;
while (inc < iter)
    [N,dNdxi]=lagrange_basis('Q4',coord);   % compute shape functions

    x = N'*nodes(:,1);
    y = N'*nodes(:,2);
    
    df1dr = dNdxi(:,1)' * nodes(:,1);
    df1ds = dNdxi(:,2)' * nodes(:,1);
    df2dr = dNdxi(:,1)' * nodes(:,2);
    df2ds = dNdxi(:,2)' * nodes(:,2);
 
    f1 = x - gpoint(1);
    f2 = y - gpoint(2);

    detF = df1dr*df2ds - df1ds*df2dr ;

    invf(1,1) =  1.0/detF * df2ds;
    invf(1,2) = -1.0/detF * df1ds;
    invf(2,1) = -1.0/detF * df2dr;
    invf(2,2) =  1.0/detF * df1dr;

    ksi = ksi - invf(1,1)*f1 - invf(1,2)*f2;
    eta = eta - invf(2,1)*f1 - invf(2,2)*f2;

    coord(1) = ksi;
    coord(2) = eta;

    if( (abs(ksi - coord(1)) < epsilon) && ...
        (abs(eta - coord(2)) < epsilon) )
        inc  = iter + 1;
        lpoint = coord;
    else
        inc = inc + 1;
    end
end