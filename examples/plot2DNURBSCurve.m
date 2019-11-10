function [sCurve,iKnot] = plot2DNURBSCurve (knots, controlPts, p, weights, noPts)

addpath ../nurbs-geopdes/inst/
addpath ../C_files/

uniqueKnots = unique(knots);
xiArr       = linspace(0,max(knots),noPts);
sCurve      = zeros(2,noPts);
iKnot       = zeros(2,length(uniqueKnots));

for i=1:noPts
  xi          = xiArr(i);   
  sCurve(1,i) = NURBSinterpolation(xi, p, knots, controlPts(:,1), weights);
  sCurve(2,i) = NURBSinterpolation(xi, p, knots, controlPts(:,2), weights);  
end

for i=1:length(uniqueKnots)
  iKnot(1,i) = NURBSinterpolation(uniqueKnots(i), p, knots, ...
                       controlPts(:,1), weights);
  iKnot(2,i) = NURBSinterpolation(uniqueKnots(i), p, knots, ...
                       controlPts(:,2), weights);
end
