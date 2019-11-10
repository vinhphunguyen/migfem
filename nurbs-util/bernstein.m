function B=bernstein(p,a,xi)

if p==0 && a==1
    B=1;
elseif p==0 && a~=1
    B=0;
else
    if a<1 || a>p+1
        B=0;
    else
        B1=bernstein(p-1,a,xi); 
        B2=bernstein(p-1,a-1,xi); 
        B=0.5*(1-xi)*B1+0.5*(1+xi)*B2;
    end
end