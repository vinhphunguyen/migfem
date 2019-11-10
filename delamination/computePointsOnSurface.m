function samPts = computePointsOnSurface(uknot,vknot,p,q,n,m,cpoints,xi,eta)

noPts  = length(xi);
samPts = zeros(noPts*noPts,3);
io = 1;
for i=1:length(eta)
    ett = eta(i);
    for j=1:length(xi)
        xii = xi(j);
        S = SurfacePoint(n,p,uknot,m,q,vknot,cpoints,3,xii,ett);
        samPts(io,:) = S;
        io = io + 1;
    end
end

