function energy = getEnergy(offsetPts,samPts)

energy = 0;

for i=1:length(offsetPts)
   p1 = offsetPts(i,:);
   p2 = samPts   (i,:);
   d  = distance(p1,p2);
   energy = energy + d*d;
end

energy = energy*0.5;