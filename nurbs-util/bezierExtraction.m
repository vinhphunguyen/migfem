function [C nb] = bezierExtraction(knot,p)
% Bezier extraction
% Taken from the master thesis of 
% Norway
% Modified by VP Nguyen in case of multiplicity = p+1

m  = length(knot)-p-1;
a  = p+1;
b  = a+1;
nb = 1;
C(:,:,1) = eye(p+1);

while b <= m
    C(:,:,nb+1) = eye(p+1);
    i=b;
    while b <= m && knot(b+1) == knot(b)
        b=b+1;
    end
    
    multiplicity = b-i+1;
    if multiplicity < p
        numerator=knot(b)-knot(a);
        for j=p:-1:multiplicity+1
            alphas(j-multiplicity)=numerator/(knot(a+j)-knot(a));
        end
        r=p-multiplicity;
        for j=1:r
            save = r-j+1;
            s = multiplicity + j;
            for k=p+1:-1:s+1
                alpha=alphas(k-s);
                C(:,k,nb)=alpha*C(:,k,nb)+(1-alpha)*C(:,k-1,nb);
            end
            if b <= m
                C(save:save+j,save,nb+1)=C(p-j+1:p+1,p+1,nb);
            end
        end
        nb=nb+1;
        if b <= m
            a=b;
            b=b+1;
        end
    elseif (multiplicity==p) || (multiplicity==p+1) % modified by VP Nguyen 
        if b <= m                                   % for case m=p+1 (delamination analysis) 
            nb=nb+1; a=b; b=b+1;
        end
    end
end

